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
# from libgl.representations import DotsRepresentation
# from libgl.representations import LinesRepresentation
# from libgl.representations import SticksRepresentation
# from libgl.representations import SpheresRepresentation
# from libgl.representations import NonBondedRepresentation
from core.vismol_config import VismolConfig
from core.vismol_selections import VismolPickingSelection as VMPick
from core.vismol_selections import VismolViewingSelection as VMSele



class VismolSession():
    """ Class doc """
    
    def __init__(self, widget=False, toolkit="gtk3", main_session=None):
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
        self.show_hide_indexes = None
        self._picking_selection_mode = False # True/False  - interchange between viewing  and picking mode
        self.selections = {"sel01": VMSele(self)}
        self.current_selection = "sel01"
        self.picking_selections = VMPick(self)
        
        if widget:
            if toolkit == "gtk3":
                from gui.vismol_gtkwidget import VismolGTKWidget
                self.selection_box_frame = None
                self.vm_widget = VismolGTKWidget(self)
                self.vm_glcore = self.vm_widget.vm_glcore
                self.vm_glcore.queue_draw()
                self.gtk_widgets_update_list = []
            if toolkit == "qt4":
                self.vm_widget = None
                raise NotImplementedError("Not implemented yet for Qt4 :(")
                quit()
        else:
            self.vm_widget = None
            raise RuntimeError("Widget not found, quitting.")
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
            self.show_hide_indexes = []
            if rep_type == "lines":
                for bond in vm_object.bonds:
                    if bond.atom_i.lines and bond.atom_j.lines:
                        self.show_hide_indexes.append(bond.atom_index_i)
                        self.show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "sticks":
                for bond in vm_object.bonds:
                    if bond.atom_i.sticks and bond.atom_j.sticks:
                        self.show_hide_indexes.append(bond.atom_index_i)
                        self.show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "ribbon":
                raise NotImplementedError("Not implementer for 'ribbon' yet.")
            
            elif rep_type == "dots":
                for atom in vm_object.atoms:
                    if atom.dots:
                        self.show_hide_indexes.append(vm_object.atoms.index(atom))
            
            elif rep_type == "nonbonded":
                for atom in vm_object.atoms:
                    if atom.nonbonded:
                        self.show_hide_indexes.append(vm_object.atoms.index(atom))
            
            elif rep_type == "impostor":
                raise NotImplementedError("Not implementer for 'impostor' yet.")
            
            elif rep_type == "spheres":
                raise NotImplementedError("Not implementer for 'spheres' yet.")
            
            elif rep_type == "surface":
                raise NotImplementedError("Not implementer for 'surface' yet.")
            
            elif rep_type == "cartoon":
                raise NotImplementedError("Not implementer for 'cartoon' yet.")
            
            if len(self.show_hide_indexes) > 0:
                vm_object.representations[rep_type].define_new_indexes_to_vbo(indexes_bonds)
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
    
    def selection_function(self, pickedID):
        """ Function doc """
        #print("selection_function")
        if pickedID is None:
            selected = None
        else:
            selected = self.atom_dic_id[pickedID]
        
        #"""     P I C K I N G     S E L E C T I O N S     """
        if self._picking_selection_mode:
            self.picking_selections.selection_function_picking(selected)
        
        else:
            self.selections[self.current_selection].selection_function_viewing(selected)
    
    def _selection_function(self, selected, _type=None, disable=True):
        #"""     P I C K I N G     S E L E C T I O N S     """
        #print("_selection_function")
        if self._picking_selection_mode:
            self.picking_selections.selection_function_picking(selected)
        
        #"""     V I E W I N G     S E L E C T I O N S     """
        else:
            self.selections[self.current_selection].selection_function_viewing(selected, _type, disable)
    




'''
    def gtk_widgets_update(self):
        """ Function doc """
        for widget in self.gtk_widgets_update_list:
            widget.update()
    
    def insert_glmenu(self, bg_menu=None, sele_menu=None, obj_menu=None, pick_menu=None):
        """ Function doc """
        
        def _viewing_selection_mode_atom(_):
            """ Function doc """
            self.viewing_selection_mode(sel_type="atom")
        
        def _viewing_selection_mode_residue(_):
            """ Function doc """
            self.viewing_selection_mode(sel_type="residue")
        
        def _viewing_selection_mode_chain(_):
            """ Function doc """
            self.viewing_selection_mode(sel_type="chain")
        
        def _selection_type_picking(_):
            """ Function doc """
            if self.selection_box_frame:
                self.selection_box_frame.change_toggle_button_selecting_mode_status(True)
            else:
                self._picking_selection_mode = True
            self.vm_glcore.queue_draw()
        
        def _selection_type_viewing(_):
            if self.selection_box_frame:
                self.selection_box_frame.change_toggle_button_selecting_mode_status(False)
            else:
                self._picking_selection_mode = False
            self.vm_glcore.queue_draw()
        
        if sele_menu is None:
            """ Standard Sele Menu """
            
            def select_test(_):
                """ Function doc """
                self.select(indexes="all")
            
            def menu_show_lines(_):
                """ Function doc """
                self.show_or_hide(rep_type="lines", show=True)
            
            def menu_hide_lines(_):
                """ Function doc """
                self.show_or_hide(rep_type="lines", show=False)
            
            def menu_show_sticks(_):
                """ Function doc """
                self.show_or_hide(rep_type="sticks", show=True)
            
            def menu_show_nonbonded(_):
                """ Function doc """
                self.show_or_hide(rep_type="nonbonded", show=True)
            
            def menu_hide_nonbonded(_):
                """ Function doc """
                self.show_or_hide(rep_type="nonbonded", show=False)
            
            def menu_hide_sticks(_):
                """ Function doc """
                self.show_or_hide(rep_type="sticks", show=False)
            
            def menu_show_spheres(_):
                """ Function doc """
                self.show_or_hide(rep_type="spheres", show=True)
            
            def menu_hide_spheres(_):
                """ Function doc """
                self.show_or_hide(rep_type="spheres", show=False)
            
            def menu_show_dots(_):
                """ Function doc """
                self.show_or_hide(rep_type="dots", show=True)
            
            def menu_hide_dots(_):
                """ Function doc """
                self.show_or_hide(rep_type="dots", show=False)
            
            def invert_selection(_):
                """ Function doc """
                self.selections[self.current_selection].invert_selection()
            
            sele_menu = { 
                    "header" : ["MenuItem", None],
                    "separator1":["separator", None],
                    "show"   : [
                                "submenu" ,{
                                            
                                            "lines"         : ["MenuItem", menu_show_lines],
                                            "sticks"        : ["MenuItem", menu_show_sticks],
                                            "spheres"       : ["MenuItem", menu_show_spheres],
                                            "dots"          : ["MenuItem", menu_show_dots],
                                            "separator2"    : ["separator", None],
                                            "nonbonded"     : ["MenuItem", menu_show_nonbonded],
                    
                                           }
                               ],
                    
                    
                    "hide"   : [
                                "submenu",  {
                                            "lines"    : ["MenuItem", menu_hide_lines],
                                            "sticks"   : ["MenuItem", menu_hide_sticks],
                                            "spheres"  : ["MenuItem", menu_hide_spheres],
                                            "dots"     : ["MenuItem", menu_hide_dots],
                                            "separator2"    : ["separator", None],
                                            "nonbonded": ["MenuItem", menu_hide_nonbonded],
                                            }
                                ],
                    
                    "Invert Selection":["MenuItem", invert_selection],
                    
                    "separator2":["separator", None],
            
                    
                    
                    "Selection type"   : [
                                "submenu" ,{
                                            
                                            "viewing"   :  ["MenuItem", _selection_type_viewing],
                                            "picking"   :  ["MenuItem", _selection_type_picking],
                                            #"separator2":["separator", None],
                                            #"nonbonded" : ["MenuItem", None],
                    
                                           }
                                        ],
                    
                    "Selection Mode"   : [
                                "submenu" ,{
                                            
                                            "Atoms"     :  ["MenuItem", _viewing_selection_mode_atom],
                                            "Residue"   :  ["MenuItem", _viewing_selection_mode_residue],
                                            "Chain"     :  ["MenuItem", _viewing_selection_mode_chain],
                                            #"separator2":["separator", None],
                                            #"nonbonded" : ["MenuItem", None],
                    
                                           }
                               ],
                    
                    "separator3":["separator", None],
                    
                    "Label Mode":  ["submenu" , {
                                            "Atom"         : [
                                                               "submenu", {
                                                                           "lines"    : ["MenuItem", None],
                                                                           "sticks"   : ["MenuItem", None],
                                                                           "spheres"  : ["MenuItem", None],
                                                                           "nonbonded": ["MenuItem", None],
                                                                           }
                                                              ],
                                            
                                            "Atom index"   : ["MenuItem", None],
                                            "residue name" : ["MenuItem", None],
                                            "residue_index": ["MenuItem", None],
                                           },
                               ]
                    }
        
        if bg_menu is None:
            """ Standard Bg Menu"""
            def open_structure_data(_):
                """ Function doc """
                self.filechooser = FileChooser()
                filename = self.filechooser.open()
                self.load(filename)
            bg_menu = { 
                    "separator0"   :["separator", None],

                    "Open File"    : ["MenuItem", open_structure_data],
                    
                    "select" : ["MenuItem", select_test],

                    "separator1":["separator", None],


                    "Selection type"   : [
                                "submenu" ,{
                                            
                                            "viewing"   :  ["MenuItem", _selection_type_viewing],
                                            "picking"   :  ["MenuItem", _selection_type_picking],
                                            #"separator2":["separator", None],
                                            #"nonbonded" : ["MenuItem", None],
                    
                                           }
                                        ],
                    
                    "Selection Mode"   : [
                                "submenu" ,{
                                            
                                            "atoms"     :  ["MenuItem", _viewing_selection_mode_atom],
                                            "residue"   :  ["MenuItem", _viewing_selection_mode_residue],
                                            "chain"     :  ["MenuItem", _viewing_selection_mode_chain],
                                            #"separator2":["separator", None],
                                            #"nonbonded" : ["MenuItem", None],
                    
                                           }
                               ],
                    
                    
                    "hide"   : [
                                "submenu",  {
                                            "lines"    : ["MenuItem", menu_hide_lines],
                                            "sticks"   : ["MenuItem", menu_hide_sticks],
                                            "spheres"  : ["MenuItem", menu_hide_spheres],
                                            "nonbonded": ["MenuItem", None],
                                            }
                                ],
                    
                    
                    "separator2":["separator", None],

                    
                    
                    "label":  ["submenu" , {
                                            "Atom"         : [
                                                               "submenu", {
                                                                           "lines"    : ["MenuItem", None],
                                                                           "sticks"   : ["MenuItem", None],
                                                                           "spheres"  : ["MenuItem", None],
                                                                           "nonbonded": ["MenuItem", None],
                                                                           }
                                                              ],
                                            
                                            "Atom index"   : ["MenuItem", None],
                                            "residue name" : ["MenuItem", None],
                                            "residue_index": ["MenuItem", None],
                                           },
                               ]
                    }

        if obj_menu is None:
            """ Standard Obj Menu"""
            obj_menu = { 
                    "OBJ menu" : ["MenuItem", None],
                    
                    
                    "separator1":["separator", None],
                    
                    
                    "show"   : [
                                "submenu" ,{
                                            
                                            "lines"    : ["MenuItem", menu_show_lines],
                                            "sticks"   : ["MenuItem", menu_show_sticks],
                                            "spheres"  : ["MenuItem", menu_show_spheres],
                                            "separator2":["separator", None],
                                            "nonbonded": ["MenuItem", None],
                    
                                           }
                               ],
                    
                    
                    "hide"   : [
                                "submenu",  {
                                            "lines"    : ["MenuItem", menu_hide_lines],
                                            "sticks"   : ["MenuItem", menu_hide_sticks],
                                            "spheres"  : ["MenuItem", menu_hide_spheres],
                                            "nonbonded": ["MenuItem", None],
                                            }
                                ],
                    
                    
                    "separator2":["separator", None],

                    
                    
                    "label":  ["submenu" , {
                                            "Atom"         : [
                                                               "submenu", {
                                                                           "lines"    : ["MenuItem", None],
                                                                           "sticks"   : ["MenuItem", None],
                                                                           "spheres"  : ["MenuItem", None],
                                                                           "nonbonded": ["MenuItem", None],
                                                                           }
                                                              ],
                                            
                                            "atomic index" : ["MenuItem", None],
                                            "residue name" : ["MenuItem", None],
                                            "residue_index": ["MenuItem", None],
                                           },
                               ]
                    }



        if pick_menu is None:
            """ Standard Sele Menu """
            pick_menu = { 
                    "header" : ["MenuItem", None],
                    
                    
                    
                    "separator1":["separator", None],
                    
                    
                    "show"   : [
                                "submenu" ,{
                                            
                                            "lines"         : ["MenuItem", menu_show_lines],
                                            "sticks"        : ["MenuItem", menu_show_sticks],
                                            "spheres"       : ["MenuItem", menu_show_spheres],
                                            "separator2"    : ["separator", None],
                                            "nonbonded"     : ["MenuItem", None],
                    
                                           }
                               ],
                    
                    
                    "hide"   : [
                                "submenu",  {
                                            "lines"    : ["MenuItem", menu_hide_lines],
                                            "sticks"   : ["MenuItem", menu_hide_sticks],
                                            "spheres"  : ["MenuItem", menu_hide_spheres],
                                            "nonbonded": ["MenuItem", None],
                                            }
                                ],
                    
                    
                    "separator2":["separator", None],

                    }
        self.vm_widget.build_glmenu(bg_menu=bg_menu, sele_menu=sele_menu,
                                    obj_menu=obj_menu, pick_menu=pick_menu)
    
    def command_line (self, entry = None):
        """ Function doc """
        cmd = entry.split()
        print (cmd)
        
        obj     = int(cmd[1]            )
        _indexes = cmd[2].split("+")
        indexes = []
        
        
        for index in _indexes:
            indexes.append(int(index))
        
        if cmd[0] == "show":
            self._show_lines (vismol_object = self.vm_objects[obj], 
                                       indexes = indexes)       
        
        if cmd[0] == "hide":
            self._hide_lines (vismol_object = self.vm_objects[obj], 
                                       indexes = indexes)  
        
        self.ctrl = True
        
        
        print (entry)
       
'''