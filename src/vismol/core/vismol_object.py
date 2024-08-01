#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import time
import numpy as np
from typing import List
from logging import getLogger
from vismol.model.atom import Atom
from vismol.model.bond import Bond
from vismol.model.chain import Chain
from vismol.model.residue import Residue
from vismol.model.molecule import Molecule
from vismol.libgl.vismol_font import VismolFont
from vismol.libgl.representations import CartoonRepresentation
from vismol.libgl.representations import DashedLinesRepresentation
from vismol.libgl.representations import DotsRepresentation
from vismol.libgl.representations import ImpostorRepresentation
from vismol.libgl.representations import LabelRepresentation
from vismol.libgl.representations import LinesRepresentation
from vismol.libgl.representations import NonBondedRepresentation
from vismol.libgl.representations import PickingDotsRepresentation
from vismol.libgl.representations import SpheresRepresentation
from vismol.libgl.representations import SticksRepresentation
from vismol.libgl.representations import SurfaceRepresentation


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
        dynamic_bonds (List): 
        c_alpha_bonds (List): 
        c_alpha_bonds (List): 
    """
    def __init__(self, vm_session: "VismolSession", index: int, name: str="UNK",
                 active: bool=False, color_palette: int=0):
        """
        Args:
            vm_session (Vismol Session): Necessary to build the
                    "atomtree_structure". The vismol_session contains the atom id
                    counter (self.vm_session.atom_id_counter).
            index (int): Unique index for the VismolObject to find it in
                    self.vm_session.vismol_objects_dic.
            name (str): Label that describes the object (default: "UNK").
            active (bool): Flag to enable/disable the object (default: False).
            color_palette (int): Number to access the color pick for carbon
                    atoms (0 = green, 1 = purple, ...) (optional).
        """
        self.vm_session = vm_session
        self.index = index
        self.name = name
        self.active = active
        self.color_palette = self.vm_session.periodic_table.get_color_palette()
        
        self.vm_config = vm_session.vm_config
        self.moving = False
        self.vm_font = VismolFont()
        self.molecule = None
        self.selected_atom_ids = set()
        self.representations = {}
        for rep_type in self.vm_config.representations_available:
            self.representations[rep_type] = None
        self.representations["ribbon_sphere"] = None
        self.representations["stick_spheres"] = None
        self.representations["labels"] = None
        self.model_mat = np.identity(4, dtype=np.float32)
        self.core_representations = {"picking_dots":None, "picking_text":None}
        self.selection_dots_vao = None
        self.selection_dot_buffer = None
        self.picking_dots_vao = None
        self.picking_dot_buffer = None
        self.colors = None
        self.color_indexes = None
        self.colors_rainbow = None
        
        self.dynamic_bonds  = [] # Pair of atoms, something like: [[0,1,1,2,3,4] , [0,1,1,2], ...]
                                # Like self.index_bonds but for each frame
        self.c_alpha_bonds = [] # List of pair of atoms defining dynamic bonds for each frame
        self.c_alpha_atoms = [] # List of pair of atoms defining C-alpha bonds
    
    
    def build_core_representations(self) -> None:
        """
        Builds core representations for the VismolObject. These representations
        are for the selection dots and the labels.
        
        """
        self.core_representations["picking_dots"] = PickingDotsRepresentation(self,
                self.vm_session.vm_glcore, active=True,
                indexes=list(self.molecule.atoms.keys()))
        self.core_representations["dash"] = DashedLinesRepresentation(self,
                self.vm_session.vm_glcore, active=True,
                indexes=self.molecule.index_bonds)
    
    
    def create_representation(self, indexes: List=None, rep_type: str="lines") -> None:
        """
        Note: The function creates a new visualization representation object
        for the VismolObject and updates its internal representations dictionary.
        It initializes the representation according to the provided rep_type
        and indexes List.
        
        Args:
            rep_type (str): The type of representation to create.
            indexes (List): List of atom indexes or bond indexes based on 
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
                    indexes=self.molecule.index_bonds)
        elif rep_type == "nonbonded":
            self.representations["nonbonded"] = NonBondedRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.non_bonded_atoms)
        elif rep_type == "impostor":
            self.representations["impostor"] = ImpostorRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()))
        elif rep_type == "sticks":
            self.representations["sticks"] = SticksRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.index_bonds)
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
                    indexes=self.molecule.index_bonds)
        elif rep_type == "ribbons":
            self.representations["ribbons"] = SticksRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.index_bonds, name="ribbons")
        elif rep_type == "ribbon_sphere":
            self.representations["ribbon_sphere"] = SpheresRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=list(self.molecule.atoms.keys()), mode=2)
        elif rep_type == "dynamic":
            self.representations["dynamic"] = SticksRepresentation(self,
                    self.vm_session.vm_glcore, active=True,
                    indexes=self.molecule.index_bonds, is_dynamic=True)
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
        _mol_atoms = self.molecule.atoms
        atom_qtty = len(_mol_atoms)
        half = int(atom_qtty/2)
        quarter = int(atom_qtty/4)
        color_step = 1.0/(atom_qtty/4.0)
        red = 0.0
        green = 0.0
        blue = 1.0
        self.colors = np.empty([atom_qtty, 3], dtype=np.float32)
        self.color_indexes = np.empty([atom_qtty, 3], dtype=np.float32)
        self.colors_rainbow = np.empty([atom_qtty, 3], dtype=np.float32)
        for i, atom in _mol_atoms.items():
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
    
    
    def set_model_matrix(self, mat: np.array) -> None:
        """
        Updates the model matrix with a copy of the numpy array provided.
        
        """
        self.model_mat = np.copy(mat)
    
