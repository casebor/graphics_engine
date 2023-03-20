#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_session.py
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

import os
import numpy as np
from vismol.utils import parser
from logging import getLogger
from vismol.core.vismol_config import VismolConfig
from vismol.core.vismol_selections import VismolPickingSelection as VMPick
from vismol.core.vismol_selections import VismolViewingSelection as VMSele

logger = getLogger(__name__)


class VismolSession():
    """ Class doc """
    
    def __init__(self, toolkit, widget=None, main_session=None):
        """ Class initialiser """
        self.main_session = None
        self.toolkit = toolkit
        self.frame = 0
        self.vm_config = VismolConfig(self)
        self.vm_objects_dic = {}
        self.atom_dic_id = {}
        self.atom_id_counter = np.uint32(0) # This variable could be replaced by len(self.atom_dic_id)?
        self.picking_selection_mode = False # True/False  - interchange between viewing  and picking mode
        self.selections = {"sel_00": VMSele(self)}
        self.current_selection = "sel_00"
        self.picking_selections = VMPick(self)
        
        self.vm_geometric_object_dic = {"pk1pk2":None, "pk2pk3":None, "pk3pk4":None}
        
        if toolkit == "Gtk_3.0":
            from vismol.gui.vismol_gtkwidget import VismolGTKWidget
            self.selection_box_frame = None
            if widget is None:
                self.vm_widget = VismolGTKWidget(self)
            else:
                self.vm_widget = widget
            self.vm_glcore = self.vm_widget.vm_glcore
            self.vm_glcore.queue_draw()
        
        elif toolkit == "Qt5":
            self.vm_widget = None
            logger.error("Not implemented yet for Qt5 :(")
            raise NotImplementedError("Not implemented yet for Qt5 :(")
            quit()
        
        else:
            self.vm_widget = None
            logger.error("Toolkit not defined or syntax error, try 'Gtk_3.0'. Quitting.")
            raise RuntimeError("Toolkit not defined or syntax error, try 'Gtk_3.0'. Quitting.")
            quit()
    
    def _add_vismol_object(self, vismol_object, show_molecule=True, autocenter=True):
        """ Function doc """
        if vismol_object.index in self.vm_objects_dic.keys():
            logger.warning("The VismolObject with id {} already exists. \
                The data will be overwritten.".format(vismol_object.index))
        self.vm_objects_dic[len(self.vm_objects_dic)] = vismol_object
        self.atom_id_counter += len(vismol_object.atoms)
        for atom in vismol_object.atoms.values():
            self.atom_dic_id[atom.unique_id] = atom
        if show_molecule:
            vismol_object.create_representation(rep_type="lines")
            vismol_object.create_representation(rep_type="nonbonded")
            if autocenter:
                self.vm_glcore.center_on_coordinates(vismol_object, vismol_object.mass_center)
            else:
                self.vm_glcore.queue_draw()
    
    def load_molecule(self, infile):
        """ Probably would be better to join this with _add_vismol_object
        """
        vismol_object, show_molecule = parser.parse_file(self, infile)
        vismol_object._generate_color_vectors(self.atom_id_counter)
        vismol_object.active = True
        self._add_vismol_object(vismol_object, show_molecule=show_molecule)
    
    def _change_attributes_for_atoms(self, atoms, rep_type, show):
        """ Function doc """
        for atom in atoms:
            try:
                if show:
                    setattr(atom, rep_type, True)
                else:
                    setattr(atom, rep_type, False)
            except AttributeError as ae:
                logger.error("Representation of type {} not implemented".format(rep_type))
                logger.error(ae)
    
    def show_or_hide(self, rep_type="lines", selection=None, show=True):
        """ Function doc """
        if selection is None:
            selection = self.selections[self.current_selection]
        
        self._change_attributes_for_atoms(selection.selected_atoms, rep_type, show)
        for vm_object in selection.selected_objects:
            if vm_object.representations[rep_type] is None:
                vm_object.create_representation(rep_type=rep_type)
            show_hide_indexes = []
            if rep_type == "lines":
                for bond in vm_object.bonds:
                    if bond.atom_i.lines and bond.atom_j.lines:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "sticks":
                for bond in vm_object.bonds:
                    if bond.atom_i.sticks and bond.atom_j.sticks:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dash":
                for bond in vm_object.bonds:
                    if bond.atom_i.dash and bond.atom_j.dash:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dynamic":
                self.define_dynamic_bonds()
                for bond in vm_object.bonds:
                    if bond.atom_i.dynamic and bond.atom_j.dynamic:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dots":
                for atom in vm_object.atoms.values():
                    if atom.dots:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "nonbonded":
                for atom in vm_object.atoms.values():
                    if atom.nonbonded:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "impostor":
                for atom in vm_object.atoms.values():
                    if atom.impostor:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "spheres":
                for atom in vm_object.atoms.values():
                    if atom.spheres:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "vdw_spheres":
                for atom in vm_object.atoms.values():
                    if atom.vdw_spheres:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "ribbons": # add by bachega at 02/02/2023
                for bond in vm_object.c_alpha_bonds:
                    if bond.atom_i.ribbons and bond.atom_j.ribbons:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
                #print('show_hide_indexes',show_hide_indexes)
                #logger.error("Not implementer for 'ribbon' yet.")
                #raise NotImplementedError("Not implementer for 'ribbon' yet.")
            
            
            elif rep_type == "surface":
                logger.error("Not implementer for 'surface' yet.")
                raise NotImplementedError("Not implementer for 'surface' yet.")
            elif rep_type == "cartoon":
                logger.error("Not implementer for 'cartoon' yet.")
                raise NotImplementedError("Not implementer for 'cartoon' yet.")
            
            if len(show_hide_indexes) > 0:
                vm_object.representations[rep_type].define_new_indexes_to_vbo(show_hide_indexes)
                vm_object.representations[rep_type].active = True
                vm_object.representations[rep_type].was_rep_ind_modified = True
                vm_object.representations[rep_type].was_sel_ind_modified = True
            else:
                vm_object.representations[rep_type].active = False
        
        self.vm_widget.queue_draw()
    
    def forward_frame(self):
        """ Function doc """
        frame = self.frame + 1
        for i, vm_object in enumerate(self.vm_objects_dic.values()):
            if frame < vm_object.frames.shape[0]:
                self.frame += 1
                self.vm_glcore.updated_coords = True
                break
            else:
                pass
        else:
            self.vm_glcore.updated_coords = False
    
    def reverse_frame(self):
        """ Function doc """
        if self.frame - 1 >= 0:
            self.frame -= 1
            self.vm_glcore.updated_coords = True
        else:
            self.vm_glcore.updated_coords = False
    
    def set_frame(self, frame=0):
        """ Function doc """
        assert frame >= 0
        self.frame = np.uint32(frame)
        self.vm_widget.queue_draw()
    
    def get_frame(self):
        """ Function doc """
        return self.frame
    
    def _selection_function_set(self, selected, _type=None, disable=True):
        """ Function doc """
        if self.picking_selection_mode: # True for picking mode
            if selected:
                assert len(selected) == 1
                selected = list(selected)[0]
                self.picking_selections.selection_function_picking(selected)
            else:
                self.picking_selections.selection_function_picking(None)
        else: # False for viewing mode
            self.selections[self.current_selection].selection_function_viewing_set(selected, _type, disable)
    
    def viewing_selection_mode(self, sel_type="atom"):
        """ Function doc """
        if self.selection_box_frame:
            self.selection_box_frame.change_sel_type_in_combobox(sel_type)
        self.selections[self.current_selection].selection_mode = sel_type
    
    def define_dynamic_bonds (self):
        """ Function doc """
        selection = self.selections[self.current_selection]
        selection_dict = {}
        vobject = None
        for atom in selection.selected_atoms:
            selection_dict[atom.atom_id] = atom
            vobject = atom.vm_object
        
        
        vobject.dynamic_bonds = []
        for frame in range(len(vobject.frames)):            
            bonds = vobject.find_bonded_and_nonbonded_atoms(selection=selection_dict, frame=frame, internal = False)
            vobject.dynamic_bonds.append(bonds)
            #print(len(bonds), bonds)
        #print(vobject.dynamic_bonds)



    # def delete_vismol_object_by_index(self, index):
    #     """ Function doc
    #     """
    #     try:
    #         vm_object = self.vm_objects_dic.pop(index)
    #     except KeyError:
    #         logger.warning("VismolObject with index {} not found".format(index))
    #         return None
    #     return vm_object
    
    # def select(self, vismol_object=None, indexes=None, sele=None):
    #     """ Function doc """
    #     if vismol_object is None:
    #         vismol_object = self.vm_objects[-1]
        
    #     if sele is None:
    #         sele = self.current_selection
        
    #     if indexes == "all":
    #         self.selections[sele].selecting_by_indexes(vismol_object=vismol_object,
    #                                                    indexes=range(0, int(len(vismol_object.atoms)/2)))
    #     self.vm_widget.queue_draw()
    
    # def disable_by_index(self, index):
    #     """ When the variable "dictionary" is active, the function accesses 
    #         a vismol object through the dictionary "self.vm_objects_dic". 
    #         Each vismol object has a unique access key (int), which, in 
    #         easyhybrid, is generated in the method: add_vismol_object.
            
    #         In the vismol interface the enable_by_index/disable_by_index methods
    #         access the vismol objects by their position in the "self.vm_objects" 
    #         list (this is because when an object is deleted in the vismol 
    #         interface, the treeview"s liststore is rewritten)
    #     """
    #     try:
    #         self.vm_objects_dic[index].active = False
    #         self.vm_glcore.queue_draw()
    #     except KeyError:
    #         logger.error("VismolObject with index {} not found".format(index))
    #         return False
    #     return True
    
    # def enable_by_index(self, index):
    #     """ When the variable "dictionary" is active, the function accesses 
    #         a vismol object through the dictionary "self.vm_objects_dic". 
    #         Each vismol object has a unique access key (int), which, in 
    #         easyhybrid, is generated in the method: add_vismol_object.
            
    #         In the vismol interface the enable_by_index/disable_by_index methods
    #         access the vismol objects by their position in the "self.vm_objects" 
    #         list (this is because when an object is deleted in the vismol 
    #         interface, the treeview"s liststore is rewritten)
    #     """
    #     try:
    #         self.vm_objects_dic[index].active = True
    #         self.vm_glcore.queue_draw()
    #     except KeyError:
    #         logger.error("VismolObject with index {} not found".format(index))
    #         return False
    #     return True
    
    # def edit_by_index(self, index):
    #     """ Function doc
    #     """
    #     try:
    #         self.vm_objects_dic[index].editing = not self.vm_objects_dic[index].editing
    #         self.vm_glcore.queue_draw()
    #     except KeyError:
    #         logger.error("VismolObject with index {} not found".format(index))
    #         return False
    #     return True
    
    # def set_color_by_index(self, vismol_object, indexes=None, color=None):
    #     """ NOT SURE WHAT THIS FUNCTION DOES
    #     """
    #     if indexes is None:
    #         indexes = []
    #     if color is None:
    #         color = np.array([0.9, 0.9, 0.9], dtype=np.float32)
        
    #     for atom_index in indexes:
    #         vismol_object.atoms[atom_index].color = color
    #     vismol_object._generate_color_vectors(do_colors=True, do_colors_idx=False,
    #                                           do_colors_raindow=False, do_vdw_dot_sizes=False,
    #                                           do_cov_dot_sizes=False)
    #     self.vm_widget.queue_draw()
    #     for rep  in vismol_object.representations.keys():
    #         if vismol_object.representations[rep]:
    #             vismol_object.representations[rep]._set_colors_to_buffer()
    #             # try:
    #             #     vismol_object.representations[rep]._set_colors_to_buffer()
    #             # except:
    #             #     print('"VisMol/vModel/Representations.py, line 123, in _set_colors_to_buffer GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ctypes.ArgumentError: argument 2: <class "TypeError">: wrong type"')
    #     return True
    
