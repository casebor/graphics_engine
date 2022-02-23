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
from utils import parser
from core.vismol_config import VismolConfig
from core.vismol_selections import VismolPickingSelection as VMPick
from core.vismol_selections import VismolViewingSelection as VMSele


class VismolSession():
    """ Class doc """
    
    def __init__(self, toolkit, widget=None, main_session=None):
        """ Class initialiser """
        self.main_session = None
        self.toolkit = toolkit
        self.frame = 0
        self.vm_vbos = []
        self.vm_config = VismolConfig(self)
        self.vm_objects_dic = {} # old Vobjects dic - include molecules
        self.vm_object_counter = 0  # Each vismol object has a unique access key (int), which is generated in the method: add_vismol_object.
        self.vm_geometric_object = []
        self.vm_geometric_object_dic = {"pk1pk2":None, "pk2pk3":None, "pk3pk4":None}
        self.atom_dic_id = {}
        self.atom_id_counter = 0
        self._picking_selection_mode = False # True/False  - interchange between viewing  and picking mode
        self.selections = {"sel01": VMSele(self)}
        self.current_selection = "sel01"
        self.picking_selections = VMPick(self)
        
        if toolkit == "gtk3":
            from gui.vismol_gtkwidget import VismolGTKWidget
            self.selection_box_frame = None
            if widget is None:
                self.vm_widget = VismolGTKWidget(self)
            else:
                self.vm_widget = widget
            self.vm_glcore = self.vm_widget.vm_glcore
            self.vm_glcore.queue_draw()
            self.gtk_widgets_update_list = []
        elif toolkit == "qt4":
            self.vm_widget = None
            raise NotImplementedError("Not implemented yet for Qt4 :(")
            quit()
        else:
            self.vm_widget = None
            raise RuntimeError("Toolkit not defined, quitting.")
            quit()
    
    def add_vismol_object(self, vismol_object=None, show_molecule=True, autocenter=True):
        """ Function doc
        """
        vismol_object.index = self.vm_object_counter
        self.vm_objects_dic[self.vm_object_counter] = vismol_object
        self.vm_object_counter += 1
        if show_molecule:
            vismol_object.create_new_representation(rep_type="lines")
            vismol_object.create_new_representation(rep_type="nonbonded")
            if autocenter:
                self.vm_glcore.center_on_coordinates(vismol_object, vismol_object.mass_center)
            else:
                self.vm_glcore.queue_draw()
    
    def load_molecule(self, infile):
        """ Function doc
        """
        vismol_object, show_molecule = parser.parse_file(infile, self)
        vismol_object.active = True
        self.add_vismol_object(vismol_object, show_molecule)
    
    def change_attributes_for_selected_atoms(self, rep_type="lines", atoms=None, show=True):
        """ Function doc """
        if atoms is None:
            atoms = []
        for atom in atoms:
            try:
                if show:
                    _r = getattr(atom, rep_type)
                    _r = True
                else:
                    _r = getattr(atom, rep_type)
                    _r = False
            except AttributeError as ae:
                print("Representation of type {} not implemented".format(rep_type))
                print(ae)
    
    def show_or_hide(self, rep_type="lines", selection=None, show=True):
        """ Function doc """
        if selection is None:
            selection = self.selections[self.current_selection]
        
        self.change_attributes_for_selected_atoms(rep_type=rep_type, atoms=selection.selected_atoms,
                                                  show=show)
        for vm_object in selection.selected_objects:
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
            
            elif rep_type == "ribbon":
                raise NotImplementedError("Not implementer for 'ribbon' yet.")
            
            elif rep_type == "dots":
                for atom in vm_object.atoms:
                    if atom.dots:
                        show_hide_indexes.append(vm_object.atoms.index(atom))
            
            elif rep_type == "nonbonded":
                for atom in vm_object.atoms:
                    if atom.nonbonded:
                        show_hide_indexes.append(vm_object.atoms.index(atom))
            
            elif rep_type == "impostor":
                raise NotImplementedError("Not implementer for 'impostor' yet.")
            
            elif rep_type == "spheres":
                raise NotImplementedError("Not implementer for 'spheres' yet.")
            
            elif rep_type == "surface":
                raise NotImplementedError("Not implementer for 'surface' yet.")
            
            elif rep_type == "cartoon":
                raise NotImplementedError("Not implementer for 'cartoon' yet.")
            
            if len(show_hide_indexes) > 0:
                vm_object.representations[rep_type].define_new_indexes_to_vbo(show_hide_indexes)
                vm_object.representations[rep_type].active = True
            else:
                vm_object.representations[rep_type].active = False
        
        self.vm_widget.queue_draw()
    
    def delete_vismol_object_by_index(self, index):
        """ Function doc
        """
        try:
            vm_object = self.vm_objects_dic.pop(index)
        except KeyError:
            print("VismolObject with index {} not found".format(index))
            return None
        return vm_object
    
    def select(self, vismol_object=None, indexes=None, sele=None):
        """ Function doc """
        if vismol_object is None:
            vismol_object = self.vm_objects[-1]
        
        if sele is None:
            sele = self.current_selection
        
        if indexes == "all":
            self.selections[sele].selecting_by_indexes(vismol_object=vismol_object,
                                                       indexes=range(0, int(len(vismol_object.atoms)/2)))
        self.vm_widget.queue_draw()
    
    def disable_by_index(self, index):
        """ When the variable "dictionary" is active, the function accesses 
            a vismol object through the dictionary "self.vm_objects_dic". 
            Each vismol object has a unique access key (int), which, in 
            easyhybrid, is generated in the method: add_vismol_object.
            
            In the vismol interface the enable_by_index/disable_by_index methods
            access the vismol objects by their position in the "self.vm_objects" 
            list (this is because when an object is deleted in the vismol 
            interface, the treeview"s liststore is rewritten)
        """
        try:
            self.vm_objects_dic[index].active = False
            self.vm_glcore.queue_draw()
        except KeyError:
            print("VismolObject with index {} not found".format(index))
            return False
        return True
    
    def enable_by_index(self, index):
        """ When the variable "dictionary" is active, the function accesses 
            a vismol object through the dictionary "self.vm_objects_dic". 
            Each vismol object has a unique access key (int), which, in 
            easyhybrid, is generated in the method: add_vismol_object.
            
            In the vismol interface the enable_by_index/disable_by_index methods
            access the vismol objects by their position in the "self.vm_objects" 
            list (this is because when an object is deleted in the vismol 
            interface, the treeview"s liststore is rewritten)
        """
        try:
            self.vm_objects_dic[index].active = True
            self.vm_glcore.queue_draw()
        except KeyError:
            print("VismolObject with index {} not found".format(index))
            return False
        return True
    
    def edit_by_index(self, index):
        """ Function doc
        """
        try:
            self.vm_objects_dic[index].editing = not self.vm_objects_dic[index].editing
            self.vm_glcore.queue_draw()
        except KeyError:
            print("VismolObject with index {} not found".format(index))
            return False
        return True
    
    def set_color_by_index(self, vismol_object, indexes=None, color=None):
        """ NOT SURE WHAT THIS FUNCTION DOES
        """
        if indexes is None:
            indexes = []
        if color is None:
            color = np.array([0.9, 0.9, 0.9], dtype=np.float32)
        
        for atom_index in indexes:
            vismol_object.atoms[atom_index].color = color
        vismol_object._generate_color_vectors(do_colors=True, do_colors_idx=False,
                                              do_colors_raindow=False, do_vdw_dot_sizes=False,
                                              do_cov_dot_sizes=False)
        self.vm_widget.queue_draw()
        for rep  in vismol_object.representations.keys():
            if vismol_object.representations[rep]:
                try:
                    vismol_object.representations[rep]._set_colors_to_buffer()
                except:
                    print('"VisMol/vModel/Representations.py, line 123, in _set_colors_to_buffer GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ctypes.ArgumentError: argument 2: <class "TypeError">: wrong type"')
        return True
    
    def set_frame(self, frame=0):
        """ Function doc """
        self.frame = frame
        self.vm_widget.queue_draw()
    
    def get_frame(self):
        """ Function doc """
        return self.frame
    
    def viewing_selection_mode(self, sel_type="atom"):
        """ Function doc
        """
        if self.selection_box_frame:
            self.selection_box_frame.change_sel_type_in_combobox(sel_type)
            
        #print(sel_type)
        self.selections[self.current_selection]._selection_mode = sel_type
    
    def _selection_function(self, selected, _type=None, disable=True):
        #"""     P I C K I N G     S E L E C T I O N S     """
        #print("_selection_function")
        if self._picking_selection_mode:
            self.picking_selections.selection_function_picking(selected)
        
        #"""     V I E W I N G     S E L E C T I O N S     """
        else:
            self.selections[self.current_selection].selection_function_viewing(selected, _type, disable)
