#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import time
import numpy as np
from logging import getLogger
from vismol.utils import matrix_operations as mop
from vismol.model.molecular_properties import solvent_dictionary as SOLV_DICT
from vismol.core.vismol_object import VismolObject
from vismol.model.atom import Atom


logger = getLogger(__name__)


class VismolViewingSelection:
    """
    TODO: Add description.

    Attributes:
        vismol_session (VismolSession): The VismolSession object to which this
                object belongs.
        active (bool): Whether or not this object is active.
        selection_mode (str): The type of selection mode used. Defines how the
                object are selected in Vismol, it can be atom, residue, chain,
                molecule. TODO: Add other selection modes?
        selected_objects (Set): A set of objects selected.
        selected_atoms (Set): A set of atoms selected. TO_REVIEW
        selected_atom_ids (Set): A set of atom ids selected. TO_REVIEW
        selected_coords (np.arra): The coordinates of the geometric center for
                the objects selected.
    """
    
    def __init__(self, vismol_session: "VismolSession"):
        """
        Args:
            vismol_session (VismolSession): The VismolSession object to which
                this object belongs.
        
        """
        self.vm_session = vismol_session
        self.active = False
        self.selection_mode = "residue"
        self.selected_objects = set()
        self.selected_atoms = set()
        self.selected_atom_ids = set()
        self.selected_coords = None
    
    def _clear_selection_buffer(self) -> None:
        """
        Clears all the atoms selected in the system if this selection is
        disabled.
        
        """
        if not self.active:
            for vm_object in self.vm_session.vm_objects_dic.values():
                for atom in vm_object.molecule.atoms.values():
                    atom.selected = False
    
    
    def _build_selected_atoms_coords_and_selected_objects_from_selected_atoms(self) -> None:
        """
        Build selected atoms' coordinates and selected objects from selected atoms.
        
        Called latter on selection_function_viewing_set
        
        """
        # Iterate through all VismolObjects in the VismolSession
        for vm_object in self.vm_session.vm_objects_dic.values():
            # Check if the current frame is within the range of frames in the VismolObject
            if self.vm_session.frame >= vm_object.molecule.frames.shape[0]:
                frame = vm_object.molecule.frames.shape[0] - 1
            else:
                frame = self.vm_session.frame
            
            # Iterate through the selected atoms in the VismolPickingSelection
            for atom in self.selected_atoms:
                # Check if the selected atom belongs to the current VismolObject
                if atom.vm_object == vm_object:
                    # Add the selected atom's ID to the selected_atom_ids set of the VismolObject
                    vm_object.selected_atom_ids.add(atom.atom_id)

        # Clear the selected_objects set in the VismolPickingSelection
        self.selected_objects.clear()

        # Iterate through the selected atoms and add their corresponding VismolObjects to the selected_objects set
        for atom in self.selected_atoms:
            self.selected_objects.add(atom.vm_object)
    
    def selection_function_viewing_set(self, selected: set, disable=True):
        """ Takes a selected atom and passes it to the appropriate selection function.
        """
        if selected is None:
            self.selected_atoms.clear()
            self.active = False
            for vm_object in self.vm_session.vm_objects_dic.values():
                vm_object.selected_atom_ids.clear()
        else:
            if self.selection_mode == "atom":
                self.selecting_by_atom(selected, disable)
            
            elif self.selection_mode == "residue":
                self.selecting_by_residue(selected, disable)
            
            elif self.selection_mode == "molecule":
                self.selecting_by_molecule(selected, disable)
            
            elif self.selection_mode == "chain":
                self.selecting_by_chain(selected, disable)
            
            elif self.selection_mode == "c-alpha":
                self.selecting_by_c_alpha_atom(selected, disable)
            
            elif self.selection_mode == "solvent":
                self.selecting_by_solvent(selected, disable)
            
            else:
                pass
            self.active = True
        self._build_selected_atoms_coords_and_selected_objects_from_selected_atoms()
    
    def selecting_by_atom(self, selected_atoms, disable=True):
        """ The "disable" variable does not allow, if the selected 
            atom is already in the selected list, to be removed. 
            
            The disable variable is "False" for when we use 
            selection by area (selection box)
        """
        self._clear_selection_buffer()
        for selected_atom in selected_atoms:
            if selected_atom in self.selected_atoms:
                if disable:
                    self.selected_atoms.discard(selected_atom)
                    selected_atom.vm_object.selected_atom_ids.discard(selected_atom.atom_id)
                    selected_atom.selected = False
            else:
                self.selected_atoms.add(selected_atom)
                selected_atom.selected = True
    
    def selecting_by_residue(self, selected_atoms, disable=True):
        """ """
        self._clear_selection_buffer()
        # if the selected atoms IS in the selected list
        for selected_atom in selected_atoms:
            if selected_atom in self.selected_atoms:
                if disable:
                    for atom in selected_atom.residue.atoms.values():
                        self.selected_atoms.discard(atom)
                        atom.vm_object.selected_atom_ids.discard(atom.atom_id)
                        atom.selected = False
            # if the selected atoms is not in the selected list add atom by atom
            else:
                for atom in selected_atom.residue.atoms.values():
                    self.selected_atoms.add(atom)
                    atom.selected = True
    
    def selecting_by_c_alpha_atom(self, selected_atoms, disable=True):
        """ """
        self._clear_selection_buffer()
        # if the selected atoms IS in the selected list
        for selected_atom in selected_atoms:
            if selected_atom in self.selected_atoms:
                if disable:
                    for atom in selected_atom.residue.atoms.values():
                        if atom.name == 'CA':
                            self.selected_atoms.discard(atom)
                            atom.vm_object.selected_atom_ids.discard(atom.atom_id)
                            atom.selected = False
                        else:
                            pass
            # if the selected atoms is not in the selected list add atom by atom
            else:
                for atom in selected_atom.residue.atoms.values():
                    if atom.name == 'CA':
                        self.selected_atoms.add(atom)
                        atom.selected = True
                    else:
                        pass
                        
    def selecting_by_molecule(self, selected_atoms, disable=True):
        """ """
        self._clear_selection_buffer()
        # if the selected atoms IS in the selected list
        for selected_atom in selected_atoms:
            if selected_atom in self.selected_atoms:
                if disable:
                    for atom in selected_atom.molecule.atoms.values():
                        self.selected_atoms.discard(atom)
                        atom.vm_object.selected_atom_ids.discard(atom.atom_id)
                        atom.selected = False
            # if the selected atoms is not in the selected list add atom by atom
            else:
                for atom in selected_atom.molecule.atoms.values():
                    self.selected_atoms.add(atom)
                    atom.selected = True
    
    def selecting_by_chain(self, selected_atoms, disable=True):
        """ 
        when disable = True
        If the object selection is disabled, all atoms in the system will be set to False 
        
        The disable variable does not allow, if the selected 
        atom is already in the selected list, to be removed. 
        
        The disable variable is "False" for when we use 
        selection by area (selection box) 
        
        """
        self._clear_selection_buffer()
        visited = set()
        for selected_atom in selected_atoms:
            if selected_atom.chain in visited:
                continue
            visited.add(selected_atom.chain)
            if selected_atom in self.selected_atoms:
                if disable:
                    for residue in selected_atom.chain.residues.values():
                        for atom in residue.atoms.values():
                            self.selected_atoms.discard(atom)
                            atom.vm_object.selected_atom_ids.discard(atom.atom_id)
                            atom.selected = False
            else:
                for residue in selected_atom.chain.residues.values():
                    for atom in residue.atoms.values():
                        self.selected_atoms.add(atom)
                        atom.selected = True

    def selecting_by_solvent(self, selected_atoms, disable=True):
        """ """
        self._clear_selection_buffer()
        # if the selected atoms IS in the selected list
        for selected_atom in selected_atoms:
            if selected_atom in self.selected_atoms:
                if disable:
                    if selected_atom.residue.name in SOLV_DICT.keys():
                        for atom in selected_atom.residue.atoms.values():
                            self.selected_atoms.discard(atom)
                            atom.vm_object.selected_atom_ids.discard(atom.atom_id)
                            atom.selected = False
            # if the selected atoms is not in the selected list add atom by atom
            else:
                if selected_atom.residue.name in SOLV_DICT.keys():
                    for atom in selected_atom.residue.atoms.values():
                        self.selected_atoms.add(atom)
                        atom.selected = True

    def selecting_by_vismol_object(self, selected_atom, disable=True):
        """ when disable = True
            If the object selection is disabled, all atoms within the 
            residue will be set to False 
            
            The disable variable does not allow, if the selected 
            atom is already in the selected list, to be removed. 
            
            The disable variable is "False" for when we use 
            selection by area (selection box)
        """
        self._clear_selection_buffer()
        visited = set()
        for selected_atom in selected_atoms:
            if selected_atom.vm_object in visited:
                continue
            visited.add(selected_atom.vm_object)
            if selected_atom in self.selected_atoms:
                if disable:
                    for atom in selected_atom.vm_object.molecule.atoms.values():
                        if atom in self.selected_atoms:
                            self.selected_atoms.remove(atom)
                            atom.selected = False
            else:
                for atom in selected_atom.vm_object.molecule.atoms.values():
                    self.selected_atoms.add(atom)
                    atom.selected = True
    
    def selecting_by_indexes(self, vismol_object, indexes, clear=False):
        """ Function doc """
        if clear:
            self._clear_selection_buffer()
        for i in indexes:
            vismol_object.atoms[i].selected = True
        
        self._build_selection_buffer()
        self._build_selected_atoms_coords_and_selected_objects_from_selected_atoms()
    
    def invert_selection(self, vismol_object=None):
        """ not workign """
        if vismol_object is None:
            for vm_object in self.vm_session.vm_objects_dic.values():
                for atom in vm_object.molecule.atoms.values():
                    atom.selected = not atom.selected
        else:
            for atom in vismol_object.molecule.atoms.values():
                atom.selected = not atom.selected
        
        self._build_selection_buffer()
        self._build_selected_atoms_coords_and_selected_objects_from_selected_atoms()
    
    def _build_selection_buffer(self):
        """ Function doc """
        self.selected_atoms.clear()
        for vm_object in self.vm_session.vm_objects_dic.values():
            for atom in vm_object.molecule.atoms.values():
                if atom.selected:
                    self.selected_atoms.add(atom)
    
    def get_selection_info(self):
        """ Function doc """
        return len(self.selected_atoms)
    


class VismolPickingSelection:
    """
             Class to manage picking selections in Vismol.
    
    Class manages the picking selections, allowing users to select up to 
    four atoms and calculate distances between selected pairs of atoms. 
    The selected atoms are stored in the picking_selections_list list, 
    and their distances are displayed when at least two atoms are selected.
    
    """
    def __init__ (self, vismol_session):
        """ 
        Initialize the VismolPickingSelection instance.

        Args:
            vismol_session (VismolSession): Reference to the VismolSession instance.
        """
        # Initialize a list to store up to 4 picking selections (selected atoms)
        self.picking_selections_list = [None] * 4

        # Initialize an empty list to keep track of the selected atoms' indexes
        self.picking_selections_list_index = []

        # Store a reference to the VismolSession instance associated with this selection
        self.vm_session = vismol_session
    
        self.distance_dict = {
                              "1-2"        : None,
                              '2-3'        : None,
                              '3-4'        : None,
                              'angle1-2-3' : None,
                              }
    
        self.vobject_picking = None
    
    
    
    
        pklist = ['pk1pk2', 'pk2pk3', 'pk3pk4']
        #self.build_pk_objects(pklist)
    
        self.xyz1 = None
        self.xyz2 = None
        self.xyz3 = None
        self.xyz4 = None
        
        self.dist_pk1_pk2 = None
        self.dist_pk2_pk3 = None
        self.dist_pk3_pk4 = None
        self.mid_point_pk1_pk2 = None
        self.mid_point_pk2_pk3 = None
        self.mid_point_pk3_pk4 = None
        self.pk_scolor = {'pk1':[0,0,0.5], 'pk2':[0,0.5,0.5], 'pk3':[0.5,0,0.5], 'pk4':[0.5,0.5,0]}
    
    def selection_function_picking(self, selected):
        """Handle picking selections.
        
        Picking selections are dynamic and versatile. They change with 
        the change of coordinates (change of frames). It is also possible 
        to relocate part of the selected atoms, for example: in a 
        selection of four atoms #1 #2 #3 #4, the user can relocate 
        only selection #3 to a different atom without destroying the 
        other selections
        
        Args:
            selected: The selected atom object or None to clear the selections.
        """
        if selected is None:
            # If selected is None, clear the list of picking selections (set all to None)
            self.picking_selections_list = [None] * len(self.picking_selections_list)
        else:
            if selected not in self.picking_selections_list:
                # If selected is not in the list, find an empty slot to store the selected atom
                for i in range(len(self.picking_selections_list)):
                    if self.picking_selections_list[i] is None:
                        self.picking_selections_list[i] = selected
                        selected = None
                        break
                if selected is not None:
                    # If there are no empty slots, replace the last element with the new selected atom
                    self.picking_selections_list[len(self.picking_selections_list) - 1] = selected
            else:
                # If selected is already in the list, remove it by setting the corresponding slot to None
                for i in range(len(self.picking_selections_list)):
                    if self.picking_selections_list[i] == selected:
                        self.picking_selections_list[i] = None
        
        # Calculate and display distances between selected atoms
        self.print_pk_distances()

        atom1 = self.picking_selections_list[0]
        atom2 = self.picking_selections_list[1]
        atom3 = self.picking_selections_list[2]
        atom4 = self.picking_selections_list[3]

        #
        '''It is necessary to check which atoms are selected (being part 
        of the list of selected atoms). The solution applied is a bit 
        brute force, but it works well.'''
        if atom1 and atom2:
            self.generate_pki_pkj_representations(vobj_label="pk1pk2", atom1=atom1, atom2=atom2)
            self.vm_session.vm_geometric_object_dic["pk1pk2"].representations["dash"].active = True
            self.xyz1 = atom1.get_coords_from_frame()
            self.xyz2 = atom2.get_coords_from_frame()
            self.dist_pk1_pk2, self.mid_point_pk1_pk2 = self.get_dist ( self.xyz1, self.xyz2)
            self.vm_session.vm_geometric_object_dic["pk1pk2"].midpoint = self.mid_point_pk1_pk2
            self.vm_session.vm_geometric_object_dic["pk1pk2"].dist = self.dist_pk1_pk2

            if atom3:
                self.xyz3 = atom3.get_coords_from_frame()

                self.get_angle_pk1_pk2_pk3()
                print ("Angle: ", self.theta_deg)
                text =  "Angle: "+ str(self.theta_deg)

                if atom4:
                    self.xyz4 = atom4.get_coords_from_frame()
                    angle = mop.dihedral(self.xyz1, self.xyz2, self.xyz3, self.xyz4)
                    print ("Dihedral: ", angle*57.297)
        else:
            if self.vm_session.vm_geometric_object_dic["pk1pk2"]:
                self.vm_session.vm_geometric_object_dic["pk1pk2"].representations["dash"].active = False

        if atom2 and atom3:
            self.xyz2 = atom2.get_coords_from_frame()
            self.xyz3 = atom3.get_coords_from_frame()
            self.generate_pki_pkj_representations(vobj_label="pk2pk3", atom1 = atom2, atom2 = atom3)
            self.vm_session.vm_geometric_object_dic["pk2pk3"].representations["dash"].active = True
            
            self.dist_pk2_pk3, self.mid_point_pk2_pk3 =  self.get_dist ( self.xyz2, self.xyz3)
            self.vm_session.vm_geometric_object_dic["pk2pk3"].midpoint = self.mid_point_pk2_pk3
            self.vm_session.vm_geometric_object_dic["pk2pk3"].dist = self.dist_pk2_pk3
            print(self.xyz2,self.xyz3,self.dist_pk2_pk3, self.mid_point_pk2_pk3)
        else:
            if self.vm_session.vm_geometric_object_dic["pk2pk3"]:
                self.vm_session.vm_geometric_object_dic["pk2pk3"].representations["dash"].active = False

        if atom3 and atom4:
            self.generate_pki_pkj_representations(vobj_label="pk3pk4", atom1=atom3, atom2=atom4)
            self.vm_session.vm_geometric_object_dic["pk3pk4"].representations["dash"].active = True
            self.dist_pk3_pk4, self.mid_point_pk3_pk4 =  self.get_dist ( self.xyz3, self.xyz4)
            self.vm_session.vm_geometric_object_dic["pk3pk4"].midpoint = self.mid_point_pk3_pk4
            self.vm_session.vm_geometric_object_dic["pk3pk4"].dist = self.dist_pk3_pk4
        else:
            if self.vm_session.vm_geometric_object_dic["pk3pk4"]:
                self.vm_session.vm_geometric_object_dic["pk3pk4"].representations["dash"].active = False
        #'''
        
        n = 1
        for atom in self.picking_selections_list:
            label = "pk"+str(n)
            print(label, atom)
            
            if atom:
                
                self.generate_pk_representations(vobj_label = label, atom = atom)
            else:
                if self.vm_session.vm_geometric_object_dic[label]:
                    self.vm_session.vm_geometric_object_dic[label].representations["picking_spheres"].active = False
                    print('here')
                else:
                    pass
                pass
            n+=1
            
    def generate_pk_representations (self, vobj_label = "pk1", atom = None):
        """ Function doc """

        if atom:
            print(atom, "<===========")
            xyz1 = atom.get_coords_from_frame(0)
            coords  = np.empty([1, 2, 3], dtype=np.float32)
            x = np.float32(xyz1[0])
            y = np.float32(xyz1[1])
            z = np.float32(xyz1[2])
            coords[0,0,:] = x, y, z
            
            vdw_rad  = atom.vdw_rad 
            cov_rad  = atom.cov_rad 
            ball_rad = atom.ball_rad
            
            self.vobject_picking = VismolObject(name                  = "pk1", 
                                                index                 = -1,
                                                vismol_session        = self.vm_session)
                                                # trajectory            = coords,
                                                # bonds_pair_of_indexes = [])
            
            self.vobject_picking.set_model_matrix(self.vm_session.vm_glcore.model_mat)
            
            atom = Atom(vismol_object = self.vobject_picking, 
                         name          ='Br'  , 
                         index         = 0,
                         residue       = None , 
                         chain         = None, 
                         atom_id       = 0,
                         occupancy     = 0, 
                         bfactor       = 0 ,
                         charge        = 0 )
            
            #atom.vdw_rad  = vdw_rad *1.05 #2.3 
            #atom.cov_rad  = cov_rad *1.05 #2.3 
            #atom.ball_rad = ball_rad*1.05 #2.3 
            atom.vdw_rad  = 2.3 
            atom.cov_rad  = 2.3 
            atom.ball_rad = 2.3 
            
            #color = [1.0,1.0,1.0]
             
            color = self.pk_scolor[vobj_label]
            #color = [1,1,1]
            atom.color = np.array(color, dtype=np.float32)
            self.vobject_picking.atoms[0] = atom
            self.vobject_picking.create_representation(rep_type="picking_spheres")
            #self.vobject_picking.representations["spheres"].color2 = [0,1,1]
            self.vm_session.vm_geometric_object_dic[vobj_label] = self.vobject_picking
            self.vm_session.vm_geometric_object_dic[vobj_label].representations["picking_spheres"].define_new_indexes_to_vbo(range(1))
            self.vm_session.vm_geometric_object_dic[vobj_label].representations["picking_spheres"].was_rep_ind_modified = True
            self.vm_session.vm_geometric_object_dic[vobj_label].representations["picking_spheres"].active = True

    def generate_pki_pkj_representations(self, vobj_label="pk1pk2", atom1=None, atom2=None):
        """ Function doc """
        xyz1 = atom1.get_coords_from_frame()
        xyz2 = atom2.get_coords_from_frame()
        frame = np.array([xyz1 , xyz2], dtype=np.float32)

        self.vobject_picking = VismolObject(name="UNK", index=-1,
                                            vismol_session                 = self.vm_session,
                                            trajectory                     = [frame],
                                            bonds_pair_of_indexes          = [0,1])
        
        '''
        The representation of distances in picking selection is 
        nothing more than a vismol object containing two 'atoms' with 
        names "pk"
        '''
        self.vobject_picking.set_model_matrix(self.vm_session.vm_glcore.model_mat)
        #'''
        
        atom1 = Atom(vismol_object = self.vobject_picking, 
                    name          ='pK'  , 
                    index         = 0,
                    residue       = None , 
                    chain         = None, 
                    atom_id       = 0,
                    occupancy     = 0, 
                    bfactor       = 0 ,
                    charge        = 0 )
        
        atom2 = Atom(vismol_object = self.vobject_picking, 
                     name          ='pK'  , 
                     index         = 1,
                     residue       = None , 
                     chain         = None, 
                     atom_id       = 1,
                     occupancy     = 0, 
                     bfactor       = 0 ,
                     charge        = 0 )

        color = [1,1,1]
        atom1.color = np.array(color, dtype=np.float32)
        atom2.color = np.array(color, dtype=np.float32)


        self.vobject_picking.atoms[0] = atom1
        self.vobject_picking.atoms[1] = atom2
        self.vobject_picking.generate_color_vectors(-1)
        
        self.vobject_picking.active      = True
        self.vobject_picking.index_bonds = [0,1]
        
        self.vobject_picking.create_representation(rep_type="dash")
        #color = [1,1,1]
        color = self.vm_session.vm_config.gl_parameters["dashed_dist_lines_color"]
        self.vobject_picking.representations["dash"].color2 = [color[0],color[1],color[2]]
        self.vm_session.vm_geometric_object_dic[vobj_label] = self.vobject_picking
    
    def update_pki_pkj_rep_coordinates (self):
        """   
        This function updates the coordinates to represent picking 
        selection distances. This is because picking selections are 
        dynamic, they change with the change of frames, hence the need 
        for this function.

        It's not optimized, but it's functional.
        """
        atom1 = self.picking_selections_list[0]
        atom2 = self.picking_selections_list[1]
        atom3 = self.picking_selections_list[2]
        atom4 = self.picking_selections_list[3]
        
        '''It is necessary to check which atoms are selected (being part 
        of the list of selected atoms). The solution applied is a bit 
        brute force, but it works well.'''
        n = 1
        for atom in self.picking_selections_list:
            label = "pk"+str(n)
            if atom:
                coords  = np.empty([1, 2, 3], dtype=np.float32)
                xyz1 = atom.get_coords_from_frame()
                x = np.float32(xyz1[0])
                y = np.float32(xyz1[1])
                z = np.float32(xyz1[2])
                coords[0,0,:] = x, y, z
                self.vm_session.vm_geometric_object_dic[label].frames = coords
                self.vm_session.vm_geometric_object_dic[label].representations["picking_spheres"].was_rep_coord_modified = True
            else:
                if self.vm_session.vm_geometric_object_dic[label]:
                    self.vm_session.vm_geometric_object_dic[label].active = False
                else:
                    pass
                pass
            n+=1


        #'''        
        if atom1 and atom2:
            self.xyz1 = atom1.get_coords_from_frame()
            self.xyz2 = atom2.get_coords_from_frame()
            coords  = np.empty([1, 2, 3], dtype=np.float32)
            x = np.float32(self.xyz1[0])
            y = np.float32(self.xyz1[1])
            z = np.float32(self.xyz1[2])
            coords[0,0,:] = x, y, z
            
            x = np.float32(self.xyz2[0])
            y = np.float32(self.xyz2[1])
            z = np.float32(self.xyz2[2])
            coords[0,1,:] = x, y, z
            
            self.vm_session.vm_geometric_object_dic['pk1pk2'].frames = coords
            self.vm_session.vm_geometric_object_dic["pk1pk2"].representations["dash"].was_rep_coord_modified = True
            #self.dist_pk1_pk2 =  ((self.xyz1[0] - self.xyz2[0])**2 + (self.xyz1[1] - self.xyz2[1])**2 + (self.xyz1[2] - self.xyz2[2])**2)**0.5
            self.dist_pk1_pk2, self.mid_point_pk1_pk2 =  self.get_dist ( self.xyz1, self.xyz2)
            self.vm_session.vm_geometric_object_dic["pk1pk2"].midpoint = self.mid_point_pk1_pk2
            self.vm_session.vm_geometric_object_dic["pk1pk2"].dist = self.dist_pk1_pk2

            if atom3:
                self.xyz3 = atom3.get_coords_from_frame()
                if atom4:
                    self.xyz4 = atom4.get_coords_from_frame()

        else:
            pass

        if atom2 and atom3:
            coords  = np.empty([1, 2, 3], dtype=np.float32)
            x = np.float32(self.xyz2[0])
            y = np.float32(self.xyz2[1])
            z = np.float32(self.xyz2[2])
            coords[0,0,:] = x, y, z
            
            x = np.float32(self.xyz3[0])
            y = np.float32(self.xyz3[1])
            z = np.float32(self.xyz3[2])
            coords[0,1,:] = x, y, z
            self.vm_session.vm_geometric_object_dic['pk2pk3'].frames = coords
            self.vm_session.vm_geometric_object_dic["pk2pk3"].representations["dash"].was_rep_coord_modified = True
            self.dist_pk2_pk3, self.mid_point_pk2_pk3 =  self.get_dist(self.xyz2, self.xyz3)
            self.vm_session.vm_geometric_object_dic["pk2pk3"].midpoint = self.mid_point_pk2_pk3
            self.vm_session.vm_geometric_object_dic["pk2pk3"].dist = self.dist_pk2_pk3

        else:
            pass
        
        if atom3 and atom4:
            coords  = np.empty([1, 2, 3], dtype=np.float32)
            x = np.float32(self.xyz3[0])
            y = np.float32(self.xyz3[1])
            z = np.float32(self.xyz3[2])
            coords[0,0,:] = x, y, z
            
            x = np.float32(self.xyz4[0])
            y = np.float32(self.xyz4[1])
            z = np.float32(self.xyz4[2])
            coords[0,1,:] = x, y, z
            self.vm_session.vm_geometric_object_dic['pk3pk4'].frames = coords
            self.vm_session.vm_geometric_object_dic["pk3pk4"].representations["dash"].was_rep_coord_modified = True
            self.dist_pk3_pk4, self.mid_point_pk3_pk4 =  self.get_dist ( self.xyz3, self.xyz4)
            self.vm_session.vm_geometric_object_dic["pk3pk4"].midpoint = self.mid_point_pk3_pk4
            self.vm_session.vm_geometric_object_dic["pk3pk4"].dist = self.dist_pk3_pk4

        else:
            pass

    def get_angle_pk1_pk2_pk3 (self):
        """ Function doc """
        self.v = [ self.xyz1[0] - self.xyz2[0], self.xyz1[1] - self.xyz2[1],   self.xyz1[2] - self.xyz2[2]]
        self.w = [ self.xyz3[0] - self.xyz2[0], self.xyz3[1] - self.xyz2[1],   self.xyz3[2] - self.xyz2[2]]
        #angle = mop.angle(self.xyz1, self.xyz3)
        
        v = np.array(self.v)
        w = np.array(self.w)
        
        dot_product = np.dot(v, w)
        magnitude_v = np.linalg.norm(v)
        magnitude_w = np.linalg.norm(w)
        cosine_theta = dot_product / (magnitude_v * magnitude_w)
        self.theta     = np.arccos(cosine_theta)
        self.theta_deg = np.degrees(self.theta)
        return self.theta_deg

    def get_dist (self, xyzi = None, xyzj = None):
        """ Function doc """
        mid_point     = [((xyzi[0]+xyzj[0])/2) , (xyzi[1]+xyzj[1])/2 , (xyzi[2]+xyzj[2])/2] 
        #print(xyzi, xyzj, mid_point)
        dist_pki_pkj = ((xyzi[0]-xyzj[0])**2 + (xyzi[1]-xyzj[1])**2 + (xyzi[2]-xyzj[2])**2)**0.5
        
        return dist_pki_pkj, mid_point

    def print_pk_distances (self):
        """ Function doc """
        # Calculate and display distances between selected atoms
        
        c = 0
        text = ''
        for atom1 in self.picking_selections_list:
            for atom2 in self.picking_selections_list[c + 1:]:
                if atom1 and atom2:
                    # Get the current frame number
                    frame = self.vm_session.get_frame()

                    # Get the coordinates of the selected atoms at the current frame (frame=None means current frame)
                    coords1 = atom1.get_coords_from_frame(frame=None)
                    coords2 = atom2.get_coords_from_frame(frame=None)

                    # Calculate the Euclidean distance between the selected atoms
                    dist = np.linalg.norm(coords1 - coords2)

                    # Get the names of the selected atoms
                    name1 = atom1.name
                    name2 = atom2.name

                    # Display the distance between the selected atoms
                    text += ' {}-{}: {:7.5f}  '.format(name1, name2,dist)
                    #print("atom", name1, "atom", name2,  dist)

            c += 1
        print(text)



