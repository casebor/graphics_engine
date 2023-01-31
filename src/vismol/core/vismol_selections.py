#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_selections.py
#  
#  Copyright 2022 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import time
import numpy as np
from vismol.utils import matrix_operations as mop
from vismol.core.vismol_object import VismolObject


class VismolViewingSelection:
    """ Class doc """
    
    def __init__(self, vismol_session):
        """ Class initialiser
        """
        self.active = False
        self.selection_mode = "residue"
        self.selected_objects = set() #dic of VisMol objects (obj)
        self.selected_atoms = set() #List of atoms objects (obj)
        self.selected_atom_ids = set() #List of atoms ids
        self.selected_coords = None #coordinate (floats) x y z
        self.vm_session = vismol_session
    
    def _clear_selection_buffer(self):
        """ If the object selection is disabled,
            all atoms in the system will be set to False
        """
        if self.active:
            pass
        else:
            for vm_object in self.vm_session.vm_objects_dic.values():
                for atom in vm_object.atoms.values():
                    atom.selected = False
                    #atom.vm_object.selected_atom_ids.discard(atom.atom_id)

    def _build_selected_atoms_coords_and_selected_objects_from_selected_atoms(self):
        """ Function doc """
        coords = []
        for vm_object in self.vm_session.vm_objects_dic.values():
            if self.vm_session.frame >= vm_object.frames.shape[0]:
                frame = vm_object.frames.shape[0] - 1
            else:
                frame = self.vm_session.frame
            for atom in self.selected_atoms:
                if atom.vm_object == vm_object:
                    # coords.append(vm_object.frames[frame][atom.atom_id])
                    vm_object.selected_atom_ids.add(atom.atom_id)
        # self.selected_coords = np.array(coords, dtype=np.float32)
        self.selected_objects.clear()
        # self.selected_atom_ids.clear()
        for atom in self.selected_atoms:
            self.selected_objects.add(atom.vm_object)
            # self.selected_atom_ids.add(atom.atom_id)
    
    def selection_function_viewing_set(self, selected, _type=None, disable=True):
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
                    for atom in selected_atom.vm_object.atoms.values():
                        if atom in self.selected_atoms:
                            self.selected_atoms.remove(atom)
                            atom.selected = False
            else:
                for atom in selected_atom.vm_object.atoms.values():
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
                for atom in vm_object.atoms.values():
                    atom.selected = not atom.selected
        else:
            for atom in vismol_object.atoms.values():
                atom.selected = not atom.selected
        
        self._build_selection_buffer()
        self._build_selected_atoms_coords_and_selected_objects_from_selected_atoms()
    
    def _build_selection_buffer(self):
        """ Function doc """
        self.selected_atoms.clear()
        for vm_object in self.vm_session.vm_objects_dic.values():
            for atom in vm_object.atoms.values():
                if atom.selected:
                    self.selected_atoms.add(atom)
    
    def get_selection_info(self):
        """ Function doc """
        return len(self.selected_atoms)
    


class VismolPickingSelection:
    """ Class doc """
    
    def __init__ (self, vismol_session):
        """ Class initialiser """
        self.picking_selections_list = [None]*4
        self.picking_selections_list_index = []
        self.vm_session = vismol_session
    
    def selection_function_picking(self, selected):
        """ Function doc """
        if selected is None:
            self.picking_selections_list = [None]*len(self.picking_selections_list)
        else:
            if selected not in self.picking_selections_list:
                for i in range(len(self.picking_selections_list)):
                    if self.picking_selections_list[i] == None:
                        self.picking_selections_list[i] = selected
                        selected = None
                        break
                if selected is not None:
                    self.picking_selections_list[len(self.picking_selections_list)-1] = selected
            else:
                for i in range(len(self.picking_selections_list)):
                    if self.picking_selections_list[i] == selected:
                        self.picking_selections_list[i] = None
        
        c = 0
        for atom1 in self.picking_selections_list:
            for atom2 in self.picking_selections_list[c+1:]:
                if atom1 and atom2:
                    frame = self.vm_session.get_frame()
                    coords1 = atom1.coords(frame)
                    coords2 = atom2.coords(frame)
                    dist = np.linalg.norm(coords1 - coords2)
                    name1 = atom1.name
                    name2 = atom2.name
                    print ("atom",name1, "atom",name2,  dist)
            c += 1
        
        atom1 = self.picking_selections_list[0]
        atom2 = self.picking_selections_list[1]
        atom3 = self.picking_selections_list[2]
        atom4 = self.picking_selections_list[3]
        
        '''
        if atom1 and atom2:
            self.refresh_pk1pk2_representations(vobj_label="pk1pk2", atom1=atom1, atom2=atom2)
            self.vm_session.vm_geometric_object_dic["pk1pk2"].representations["dash"].active = True
            if atom3:
                xyz1 = atom1.coords()
                xyz2 = atom2.coords()
                xyz3 = atom3.coords()
                
                xyz1 = [ xyz1[0] - xyz2[0], xyz1[1] - xyz2[1],   xyz1[2] - xyz2[2]]
                xyz3 = [ xyz3[0] - xyz2[0], xyz3[1] - xyz2[1],   xyz3[2] - xyz2[2]]
                angle = mop.angle(xyz1, xyz3)
                
                print ("Angle: ", angle*57.297)
                text =  "Angle: "+ str( angle*57.297)
                self.vm_session.main_session.statusbar_main.push(1,text)
                if atom4:
                    xyz4 = atom4.coords()
                    angle = mop.dihedral(xyz1, xyz2, xyz3, xyz4)
                    print ("Dihedral: ", angle*57.297)
        
        else:
            if self.vm_session.vm_geometric_object_dic["pk1pk2"]:
                self.vm_session.vm_geometric_object_dic["pk1pk2"].representations["dash"].active = False
        if atom2 and atom3:
            self.refresh_pk1pk2_representations(vobj_label="pk2pk3", atom1=atom2, atom2=atom3)
            self.vm_session.vm_geometric_object_dic["pk2pk3"].representations["dash"].active = True
        else:
            if self.vm_session.vm_geometric_object_dic["pk2pk3"]:
                self.vm_session.vm_geometric_object_dic["pk2pk3"].representations["dash"].active = False
        
        if atom3 and atom4:
            self.refresh_pk1pk2_representations(vobj_label="pk3pk4", atom1=atom3, atom2=atom4)
            self.vm_session.vm_geometric_object_dic["pk3pk4"].representations["dash"].active = True
        else:
            if self.vm_session.vm_geometric_object_dic["pk3pk4"]:
                self.vm_session.vm_geometric_object_dic["pk3pk4"].representations["dash"].active = False
        '''
    def refresh_pk1pk2_representations(self, vobj_label="pk1pk2", atom1=None, atom2=None):
        """ Function doc """
        xyz1 = atom1.coords()
        xyz2 = atom2.coords()
        frame = np.array(xyz1 + xyz2, dtype=np.float32)
        
        if self.vm_session.vm_geometric_object_dic[vobj_label]:
            self.vm_session.vm_geometric_object_dic[vobj_label].frames = [frame]
            self.vm_session.vm_geometric_object_dic[vobj_label].representations["dash"]._make_gl_vao_and_vbos()
            self.vm_session.vm_geometric_object_dic[vobj_label].active = True
        else:
            atoms = []
            atoms.append({
                          "index"      : 0             , 
                          "name"       : "pK"           , 
                          "resi"       : ""            , 
                          "resn"       : ""            , 
                          "chain"      : ""            , 
                          "symbol"     : "pK"           , 
                          "occupancy"  : 00.00         , 
                          "bfactor"    : 0.00          , 
                          "charge"     : 0.00           
                          })
            atoms.append({
                          "index"      : 1     , 
                          "name"       : "pK"   , 
                          "resi"       : ""    , 
                          "resn"       : ""    , 
                          "chain"      : ""    , 
                          "symbol"     : "pK"   , 
                          "occupancy"  : 00.00 , 
                          "bfactor"    : 0.00  , 
                          "charge"     : 0.00   
                          })
            
            frame = np.array(xyz1 + xyz2, dtype=np.float32)
            self.vobject_picking = VismolObject(name="UNK", index=-1,
                                                vismol_session                 = self.vm_session,
                                                trajectory                     = [frame],
                                                bonds_pair_of_indexes          = [0,1])
            
            self.vobject_picking.active = True
            self.vobject_picking.set_model_matrix(self.vm_session.vm_glcore.model_mat)
            
            #self.vobject_picking.create_representation(rep_type = "dash")
            self.vm_session.vm_geometric_object_dic[vobj_label] = self.vobject_picking

