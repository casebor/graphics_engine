#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import time
import numpy as np
from logging import getLogger
from vismol.model.atom import Atom
from vismol.model.bond import Bond
from vismol.model.chain import Chain
from vismol.model.residue import Residue
from vismol.model.molecule import Molecule
from vismol.model.topology import Topology
from vismol.libgl.vismol_font import VismolFont
from vismol.libgl.representations import DotsRepresentation
from vismol.libgl.representations import LabelRepresentation
from vismol.libgl.representations import LinesRepresentation
from vismol.libgl.representations import SticksRepresentation
from vismol.libgl.representations import CartoonRepresentation
from vismol.libgl.representations import SpheresRepresentation
from vismol.libgl.representations import SurfaceRepresentation
from vismol.libgl.representations import ImpostorRepresentation
from vismol.libgl.representations import NonBondedRepresentation
from vismol.libgl.representations import DashedLinesRepresentation
from vismol.libgl.representations import PickingDotsRepresentation
import vismol.utils.molecular_operations as molop


logger = getLogger(__name__)


class VismolObject:
    """
    A VismolObject contains the information necessary for openGL to draw 
    a model on the screen. Everything that is represented in the graphical 
    form is stored in the form of a VismolObject.
    
    Attributes:
        vm_session (Vismol Session): Necessary to build the
                "atomtree_structure". The vismol_session contains the atom id
                counter (self.vm_session.atom_id_counter).
        index (int): Unique index for the VismolObject to find it in
                self.vm_session.vismol_objects_dic.
        name (str): Label that describes the object (default: "UNK").
        active (bool): Flag to enable/disable the object (default: False).
        color_palette (int): Number to access the color pick for carbon
                atoms (0 = green, 1 = purple, ...) (optional).
        vm_config (VismolConfig): Simple configuration object containing the
                information needed to represent the Vismol object.
        moving (bool): Flag to indicate whether or not the object is being
                modified in the window. This is not the type of editing that
                changes the structure of the Vismol object, but rather its
                matrices.
        vm_font (VismolFont): A VismolFont object used to render text.
        molecule (Molecule): A VisMol object contains one molecule. The hierarchy of
                the model is: VisMol object -> Molecule -> Chain -> Residue -> Atom
        representations (Dict): A dictionary with the representations available
                for this VisMol object. When creating the VisMol object, all of
                the representations are None.
        model_mat (np.array): A numpy array containing the matrix of the model
                used for the visualization. This matrix is modified depending
                on the rotation, translation and scale in the window.
        core_representations (Dict): A dictionary containing the representation
                for picking elements in the window.
        selection_dots_vao (OpenGL.VertexArray): OpenGL object necessary for
                rendering the selections.
        selection_dots_buffer (OpenGL.Buffer): OpenGL obejct necessary for
                rendering the selections.
        picking_dots_vao (OpenGL.VertexArray): OpenGL obejct necessary for
                rendering the text in selections.
        picking_dots_buffer (OpenGL.Buffer): OpenGL obejct necessary for
                rendering the text in selections.
        colors (np.array): Numpy array containing the color information of all
                the atoms of the molecule. This usually corresponds to the
                element color of each atom.
        color_indexes (np.array): Numpy array containing the color id for each
                atom in the Vismol session. This color is used in the OpenGL
                engine to identify which atoms are selected with the mouse. The
                color id is unique for a single atom in the whole session.
        colors_rainbow (np.array): Numpy array containing the color coding for
                atoms in the rainbow scheme.
        dynamic_bonds (list): 
        c_alpha_bonds (list): 
        c_alpha_bonds (list): 
    """
    def __init__(self, vismol_session: "VismolSession", index: int, name: str="UNK",
                 active: bool=False):
        """
        Args:
            vismol_session (Vismol Session): Necessary to build the
                    "atomtree_structure". The vismol_session contains the atom id
                    counter (self.vismol_session.atom_id_counter).
            index (int): Unique index for the VismolObject to find it in
                    self.vismol_session.vismol_objects_dic.
            name (str): Label that describes the object (default: "UNK").
            active (bool): Flag to enable/disable the object (default: False).
            color_palette (int): Number to access the color pick for carbon
                    atoms (0 = green, 1 = purple, ...) (optional).
        """
        self.vm_session = vismol_session
        self.index = index
        self.name = name
        self.active = active
        self.color_palette = self.vm_session.periodic_table.get_color_palette()
        self.vm_config = self.vm_session.vm_config
        # TODO: this attribute is not being used anywhere
        self.moving = False
        self.vm_font = VismolFont()
        self.molecule = None
        self.selected_atom_ids = set()
        self.representations = {}
        for rep_type in self.vm_config.representations_available:
            self.representations[rep_type] = None
        # TODO: Why are these three representations not included in the
        #       configuration file?
        self.representations["ribbon_sphere"] = None
        self.representations["stick_spheres"] = None
        self.representations["labels"] = None
        self.model_mat = np.identity(4, dtype=np.float32)
        # TODO: The name of this attribute is a bit confusing, better change it
        self.core_representations = {"picking_dots":None, "picking_text":None}
        # TODO: What is the difference between these colors?
        #       Improve the naming or add an explanation as comment
        self.colors = None
        self.color_indexes = None
        self.colors_rainbow = None
        # TODO: We should create the topology in the parser and here only use it,
        #       not define it. Should this be in the molecule instead of here?
        #       This can be only a dictionary of unbonded, single, double,
        #       triple bonds and aromatic. Probable better make a new class?
        #       {"unbonded": [at1, at5, at32, ...],
        #        "single": [[at2, at3], [at7, at8], ...],
        #        "double": [[at2, at4], [at7, at9], ...],
        #        "triple": [[at11, at13], [at21, at29], ...],
        #        "aromatic": [[at22, at23, at24, at25, at26, at27], ...]}
        self.topology = None
        # TODO: This can be outside topology since it is more for visualization
        self.dynamic_bonds  = [] # Pair of atoms, something like: [[0,1,1,2,3,4] , [0,1,1,2], ...]
                                # Like self.index_bonds but for each frame
        # TODO: Should these two be inside the topology?
        self.c_alpha_bonds = [] # list of pair of atoms defining dynamic bonds for each frame
        self.c_alpha_atoms = [] # list of pair of atoms defining C-alpha bonds
    
    def build_core_representations(self) -> None:
        """
        Builds core representations for the VismolObject. These representations
        are for the selection dots and the labels.
        
        """
        from vismol.gui.vismol_widget_main import VismolWidgetMain
        if isinstance(self.vm_session.vm_widget, VismolWidgetMain):
            self.core_representations["picking_dots"] = PickingDotsRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()))
            self.core_representations["dash"] = DashedLinesRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.topology.bonds_pair_of_indexes)
        else:
            logger.info("Function build_core_representations should be called "\
                "only within a VismolWidgetMain instance")
    
    
    def create_representation(self, indexes: list=None, rep_type: str="lines") -> None:
        """
        Note: The function creates a new visualization representation object
        for the VismolObject and updates its internal representations dictionary.
        It initializes the representation according to the provided rep_type
        and indexes list.
        
        Args:
            rep_type (str): The type of representation to create.
            indexes (list): list of atom indexes or bond indexes based on 
                    the representation type.
        
        Raises:
            NotImplementedError (Error): If the representation type is not
                    implemented, the program should raise an error.
        
        """
        if rep_type == "dots":
            self.representations["dots"] = DotsRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()))
        elif rep_type == "lines":
            self.representations["lines"] = LinesRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.topology.bonds_pair_of_indexes)
                    # indexes=self.molecule.index_bonds)
        elif rep_type == "nonbonded":
            self.representations["nonbonded"] = NonBondedRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.topology.non_bonded_atoms)
        elif rep_type == "impostor":
            self.representations["impostor"] = ImpostorRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()))
        elif rep_type == "sticks":
            self.representations["sticks"] = SticksRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.topology.bonds_pair_of_indexes)
        elif rep_type == "stick_spheres":
            self.representations["stick_spheres"] = SpheresRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()), mode=3)
        elif rep_type == "spheres":
            self.representations["spheres"] = SpheresRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()))
        elif rep_type == "picking_spheres":
            self.representations["picking_spheres"] = SpheresRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()), mode=1)
        elif rep_type == "vdw_spheres":
            self.representations["vdw_spheres"] = SpheresRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()), vdw=True)
        elif rep_type == "dash":
            self.representations["dash"] = DashedLinesRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.topology.bonds_pair_of_indexes)
        elif rep_type == "ribbons":
            self.representations["ribbons"] = SticksRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.topology.bonds_pair_of_indexes, name="ribbons")
        elif rep_type == "ribbon_sphere":
            self.representations["ribbon_sphere"] = SpheresRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()), mode=2)
        elif rep_type == "dynamic":
            self.representations["dynamic"] = SticksRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.topology.bonds_pair_of_indexes, is_dynamic=True)
        elif rep_type == "labels":
            self.representations["labels"] = LabelRepresentation(self,
                    self.vm_session.vm_glcore, indexes=[0,1,2], labels=None,
                    color=[1, 1, 0, 1])
        elif rep_type == "surface":
            self.representations["surface"] = SurfaceRepresentation(self,
                    self.vm_session.vm_glcore, active=True, indexes=[],
                    name="surface", is_dynamic=False)
        elif rep_type == "cartoon":
            self.representations["cartoon"] = CartoonRepresentation(self,
                    self.vm_session.vm_glcore, active=True, indexes=[],
                    name="cartoon", rep_type="mol")
        else:
            logger.error("Representation {} not implemented".format(rep_type))
            raise NotImplementedError("Representation {} not implemented".format(rep_type))
    
    def generate_color_vectors(self) -> None:
        """
        This method generates and collect the color information for the atoms
        of the Vismol object. The rainbow colors are generated, while the colors
        and the color ids are gathered from the Atom objects. The method assigns
        to each atom of the system a unique identifier. This identifier is used
        in the selection function. There are no two atoms with the same color ID
        in the whole Vismol session.
        
        """
        mol_atoms = self.molecule.atoms
        atom_qtty = len(mol_atoms)
        half = int(atom_qtty/2)
        quarter = int(atom_qtty/4)
        color_step = 1.0/(atom_qtty/4.0)
        red = 0.0
        green = 0.0
        blue = 1.0
        self.colors = np.empty([atom_qtty, 3], dtype=np.float32)
        self.color_indexes = np.empty([atom_qtty, 3], dtype=np.float32)
        self.colors_rainbow = np.empty([atom_qtty, 3], dtype=np.float32)
        for i, atom in mol_atoms.items():
            atom.generate_atom_unique_color_id()
            self.colors[i] = atom.color
            self.color_indexes[i] = atom.color_id
            if i <= 1*quarter:
                self.colors_rainbow[i,:] = red, green, blue
                green += color_step
            if (i >= 1*quarter) and (i <= 2*quarter):
                self.colors_rainbow[i,:] = red, green, blue
                blue -= color_step
            if (i >= 2*quarter) and (i <= 3*quarter):
                self.colors_rainbow[i,:] = red, green, blue
                red += color_step
            if (i >= 3*quarter) and (i <= 4*quarter):
                self.colors_rainbow[i,:] = red, green, blue
                green -= color_step
    
    def set_model_matrix(self, in_mat: np.array) -> None:
        """
        Updates the model matrix with a copy of the numpy array provided.
        
        """
        self.model_mat = np.copy(in_mat)
    
    def add_new_atom(self, atom: Atom, unique_id: int) -> None:
        """ Docs
        """
        if self.molecule is None:
            molecule = Molecule(self, name="Molecule")
            topo = Topology(molecule=molecule)
            molecule.topology = topo
            self.molecule = molecule
            chain = Chain(self, name="BX", molecule=molecule)
            molecule.chains[chain.name] = chain
            residue = Residue(self, name="RX", index=1, chain=chain)
            chain.residues[residue.index] = residue
            atom.residue = residue
            atom.chain = chain
            atom.molecule = molecule
            residue.atoms[atom.atom_id] = atom
            molecule.atoms[atom.atom_id] = atom
            molecule.frames = np.zeros([1,1,3], dtype=np.float32)
            molecule.frames[0,0] = atom.coords
            molecule.cov_radii_array = np.array([atom.cov_rad], dtype=np.float32)
            index_bonds = molop.get_bond_pairs(molecule=molecule,
                                atom_ids=np.array([atom.atom_id], dtype=np.int32),
                                cov_radii_array=np.array([atom.cov_rad], dtype=np.float32),
                                gridsize=1.2, maxbond=2.4, tolerance=1.4)
            bonds_from_pairs = molop.bonds_from_pair_of_indexes_list(molecule=molecule,
                                                    index_bonds=index_bonds)
            topo.bonds_pair_of_indexes = bonds_from_pairs
            nonbonded_atoms = molop.get_non_bonded_from_bonded_list(molecule=molecule,
                                                            index_bonds=bonds_from_pairs)
            topo.non_bonded_atoms = nonbonded_atoms
            
        else:
            residue = self.molecule.chains["BX"].residues[1]
            atom.residue = residue
            atom.chain = self.molecule.chains["BX"]
            atom.molecule = self.molecule
            residue.atoms[atom.atom_id] = atom
            self.molecule.atoms[atom.atom_id] = atom
            frame = np.vstack((self.molecule.frames[0], np.array([atom.coords])))
            self.molecule.frames = np.empty([1, frame.shape[0], 3], dtype=np.float32)
            self.molecule.frames[0] = frame
            self.molecule.cov_radii_array = np.hstack((self.molecule.cov_radii_array, atom.cov_rad))
            index_bonds = molop.get_bond_pairs(molecule=self.molecule,
                                atom_ids=np.array(list(self.molecule.atoms.keys()), dtype=np.int32),
                                cov_radii_array=self.molecule.cov_radii_array,
                                gridsize=1.2, maxbond=2.4, tolerance=1.4)
            bonds_from_pairs = molop.bonds_from_pair_of_indexes_list(molecule=self.molecule,
                                                    index_bonds=index_bonds)
            self.molecule.topology.bonds_pair_of_indexes = bonds_from_pairs
            nonbonded_atoms = molop.get_non_bonded_from_bonded_list(molecule=self.molecule,
                                                            index_bonds=bonds_from_pairs)
            self.molecule.topology.non_bonded_atoms = nonbonded_atoms
    
