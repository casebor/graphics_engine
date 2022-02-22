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
import time
import numpy as np
from pprint import pprint
from utils import parser
import libgl.shapes as shapes
from core.vismol_config import VismolConfig
from core.vismol_selections import VismolPickingSelection as VMPick
from core.vismol_selections import VismolViewingSelection as VMSele
from libgl.representations import DotsRepresentation
from libgl.representations import LinesRepresentation
from libgl.representations import SticksRepresentation
from libgl.representations import SpheresRepresentation
from libgl.representations import NonBondedRepresentation
from gui.vismol_widget import VismolWidget
from gui.gtk_widgets.player import PlayerFrame
from gui.gtk_widgets.filechooser import FileChooser
from gui.gtk_widgets.vismol_tools import VismolStatusBar
from gui.gtk_widgets.vismol_tools import VismolGoToAtomWindow2
from gui.gtk_widgets.vismol_tools import VismolTrajectoryFrame
from gui.gtk_widgets.vismol_tools import VismolSelectionTypeBox


class ShowHideVisMol:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
    
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
    
    def show_or_hide_by_object(self, rep_type="lines", vismol_object=None, selection_table=None, show=True):
        """ Function doc """
        if selection_table is None:
            selection_table = []
        
        atoms = []
        for atom_index in selection_table:
            atoms.append(vismol_object.atoms[atom_index])
        self.change_attributes_for_selected_atoms(rep_type=rep_type, atoms=atoms, show=show)
        
        if rep_type in ["lines", "sticks", "ribbon"]:
            indexes_bonds = []
            for bond in vismol_object.bonds:
                if rep_type == "lines":
                    if bond.atom_i.lines and bond.atom_j.lines:
                        indexes_bonds.append(bond.atom_index_i)
                        indexes_bonds.append(bond.atom_index_j)
                    else:
                        pass
                elif rep_type == "sticks":
                    if bond.atom_i.sticks and bond.atom_j.sticks:
                        indexes_bonds.append(bond.atom_index_i)
                        indexes_bonds.append(bond.atom_index_j)
                    else:
                        pass
                elif rep_type == "ribbon":
                    raise NotImplementedError("Not implementer for ribbon yet.")
            
            if vismol_object.representations[rep_type] is None:
                if len(indexes_bonds) > 0:
                    rep = SticksRepresentation(name=rep_type, active=True,
                                               vismol_object=vismol_object,
                                               vismol_glcore=self.vm_widget.vm_glcore,
                                               indexes=indexes_bonds)
                    vismol_object.representations[rep_type] = rep
            else:
                if len(indexes_bonds) > 0:
                    indexes_bonds = np.array(indexes_bonds, dtype=np.uint32)
                    vismol_object.representations[rep_type].define_new_indexes_to_vbo(indexes_bonds)
                    vismol_object.representations[rep_type].active = True
                else:
                    vismol_object.representations[rep_type].active = False
        else:
            indexes = []
            if rep_type == "dots":
                for atom in vismol_object.atoms:
                    if atom.dots:
                        index = vismol_object.atoms.index(atom)
                        indexes.append(index)
                
                if vismol_object.representations[rep_type] is None:
                    rep = DotsRepresentation(name=rep_type, active=True,
                                             vismol_object=vismol_object,
                                             vismol_glcore=self.vm_widget.vm_glcore,
                                             indexes=indexes)
                    vismol_object.representations[rep.name] = rep
                else:
                    if len(indexes) > 0:
                        indexes = np.array(indexes, dtype=np.uint32)
                        vismol_object.representations[rep_type].define_new_indexes_to_vbo(indexes)
                        vismol_object.representations[rep_type].active = True
                    else:
                        vismol_object.representations[rep_type].active = False
            
            if rep_type == "nonbonded":
                for atom in vismol_object.atoms:
                    if atom.nonbonded:
                        index = vismol_object.atoms.index(atom)
                        indexes.append(index)
                
                if vismol_object.representations[rep_type] is None:
                    rep = NonBondedRepresentation(name=rep_type, active=True,
                                                  vismol_object=vismol_object,
                                                  vismol_glcore=self.vm_widget.vm_glcore,
                                                  indexes=indexes)
                    vismol_object.representations[rep_type] = rep
                else:
                    if len(indexes) > 0:
                        indexes = np.array(indexes, dtype=np.uint32)
                        vismol_object.representations[rep_type].define_new_indexes_to_vbo(indexes)
                        vismol_object.representations[rep_type].active = True
                    else:
                        vismol_object.representations[rep_type].active = False
            
            if  rep_type == "spheres":
                for atom in vismol_object.atoms:
                    if atom.spheres:
                        index = vismol_object.atoms.index(atom)
                        indexes.append(index)
                
                if vismol_object.representations["spheres"] is None:
                    if len(indexes) > 0:
                        rep = SpheresRepresentation(name=rep_type, active=True,
                                                    vismol_object=vismol_object,
                                                    vismol_glcore=self.vm_widget.vm_glcore,
                                                    indexes=indexes)
                        rep._create_sphere_data()
                        vismol_object.representations[rep_type] = rep
                else:
                    if len(indexes) > 0:
                        vismol_object.representations[rep_type].update_atomic_indexes(indexes=indexes)
                    else:
                        vismol_object.representations[rep_type].active = False
        self.vm_widget.queue_draw()
    
    def show_or_hide(self, rep_type="lines", selection=None, show=True):
        """ Function doc """
        if selection is None:
            selection = self.selections[self.current_selection]
        
        self.change_attributes_for_selected_atoms(rep_type=rep_type, atoms=selection.selected_atoms,
                                                  show=show)
        for vismol_object in selection.selected_objects:
            if rep_type in ["lines","sticks","ribbon"]:
                indexes_bonds = []
                for bond in vismol_object.bonds:
                    if rep_type == "lines":
                        if bond.atom_i.lines  and  bond.atom_j.lines:
                            indexes_bonds.append(bond.atom_index_i)
                            indexes_bonds.append(bond.atom_index_j)
                    
                    if rep_type == "sticks":
                        if bond.atom_i.sticks  and  bond.atom_j.sticks:
                            indexes_bonds.append(bond.atom_index_i)
                            indexes_bonds.append(bond.atom_index_j)
                
                if vismol_object.representations[rep_type] is None:
                    if len(indexes_bonds) > 0:
                        rep = SticksRepresentation(name=rep_type, active=True,
                                                   vismol_object=vismol_object,
                                                   vismol_glcore=self.vm_widget.vm_glcore,
                                                   indexes=indexes)
                        vismol_object.representations[rep.name] = rep
                else:
                    if len(indexes_bonds) == 0:
                        vismol_object.representations[rep_type].active = False
                    else:
                        indexes_bonds = np.array(indexes_bonds, dtype=np.uint32)
                        vismol_object.representations[rep_type].define_new_indexes_to_vbo(indexes_bonds)
                        vismol_object.representations[rep_type].active = True
            
            else:
                indexes = []
                if rep_type == "dots":
                    indexes = []
                    for atom in vismol_object.atoms:
                        if atom.dots:
                            index = vismol_object.atoms.index(atom)
                            indexes.append(index)
                    if vismol_object.representations[rep_type] is None:
                        rep  = DotsRepresentation(name=rep_type, active=True,
                                                  vismol_object=vismol_object,
                                                  vismol_glcore=self.vm_widget.vm_glcore,
                                                  indexes=indexes)
                        vismol_object.representations[rep.name] = rep 
                    else:
                        if len(indexes) == 0:
                            vismol_object.representations[rep_type].active = False
                        else:
                            indexes = np.array(indexes, dtype=np.uint32)
                            vismol_object.representations[rep_type].define_new_indexes_to_vbo(indexes)
                            vismol_object.representations[rep_type].active = True
                
                if rep_type == "nonbonded":
                    indexes = []
                    for atom in vismol_object.atoms:
                        if atom.nonbonded:
                            index = vismol_object.atoms.index(atom)
                            indexes.append(index)
                    if vismol_object.representations[rep_type] is None:
                        rep  = NonBondedRepresentation(name=rep_type, active=True,
                                                       vismol_object=vismol_object,
                                                       vismol_glcore=self.vm_widget.vm_glcore,
                                                       indexes=indexes)
                        vismol_object.representations[rep.name] = rep 
                    else:
                        if len(indexes) == 0:
                            vismol_object.representations[rep_type].active = False
                        else:
                            indexes = np.array(indexes, dtype=np.uint32)
                            vismol_object.representations[rep_type].define_new_indexes_to_vbo(indexes)
                            vismol_object.representations[rep_type].active = True
                
                if  rep_type == "spheres":
                    atoms2spheres = []
                    for atom in vismol_object.atoms:
                        if atom.spheres:
                            atoms2spheres.append(atom)
                            index = vismol_object.atoms.index(atom)
                            indexes.append(index)
                    if vismol_object.representations["spheres"] is None:
                        if len(atoms2spheres) > 0:
                            rep  = SpheresRepresentation(name=rep_type, active=True,
                                                         vismol_object=vismol_object,
                                                         vismol_glcore=self.vm_widget.vm_glcore,
                                                         indexes=indexes)
                            rep._create_sphere_data()
                            vismol_object.representations[rep.name] = rep
                    else:
                        if len(atoms2spheres) == 0:
                            vismol_object.representations[rep_type].active = False
                        else:
                            vismol_object.representations[rep_type].update_atomic_indexes(indexes=indexes)
        self.vm_widget.queue_draw()

class VismolSession(ShowHideVisMol):
    """ Class doc """
    
    def __init__(self, glwidget=False, toolkit="gtk3", main_session=None):
        """ Class initialiser """
        self.main_session = None
        self.toolkit = toolkit
        self.vm_config = VismolConfig(self)
        self.vm_objects_dic = {} # old Vobjects dic - include molecules
        self.vm_object_counter = 0  # Each vismol object has a unique access key (int), which is generated in the method: add_vismol_object_to_vismol_session.
        self.vm_geometric_object = []
        self.vm_vbos = []
        self.vm_geometric_object_dic = {"pk1pk2": None, "pk2pk3": None, "pk3pk4": None}
        self.atom_id_counter = 0
        self.atom_dic_id = {}
        self._picking_selection_mode = False # True/False  - interchange between viewing  and picking mode
        self.frame = 0
        #---------------------------------------------------------------
        #  VIEWING SELECTIONS
        #---------------------------------------------------------------
        self.selections = {"sel01": VMSele(self)}
        self.current_selection = "sel01"
        #---------------------------------------------------------------
        #  PICKING SELECTIONS
        #---------------------------------------------------------------
        self.picking_selections = VMPick(self)
        #---------------------------------------------------------------------------
        # GTK WIDGETS
        #---------------------------------------------------------------------------
        if glwidget:
            if toolkit == "gtk3":
                self.selection_box_frame = None
                self.vm_widget = VismolWidget(self)
                self.vm_glcore = self.vm_widget.vm_glcore
                self.vm_glcore.queue_draw()
                self.gtk_widgets_update_list = []
            if toolkit == "qt4":
                raise NotImplementedError()
        else:
            self.vm_widget = None
    
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


    def add_vismol_object_to_vismol_session(self, vismol_object=None, show_molecule=True,
                                            autocenter=True):
        """ Function doc """
        vismol_object.index = self.vm_object_counter
        self.vm_objects_dic[self.vm_object_counter] = vismol_object
        self.vm_object_counter += 1
        # self.append_vismol_object_to_vismol_objects_list_store(vismol_object)
        
        if show_molecule:
            vismol_object.create_new_representation(rep_type="lines")
            vismol_object.create_new_representation(rep_type="nonbonded")
            if autocenter:
                self.vm_glcore.center_on_coordinates(vismol_object, vismol_object.mass_center)
            else:
                self.vm_glcore.queue_draw()
            self.gtk_widgets_update()
    
    def load(self, infile):
        """ Function doc
        """
        vismol_object, show_molecule = parser.parse_file(infile, self)
        vismol_object.active = True
        self.add_vismol_object_to_vismol_session(vismol_object, show_molecule)
    
    # def append_vismol_object_to_vismol_objects_list_store(self, vismol_object):
    #     """ This function adds new structures to "vm_objects_list_store". 
    #         The vm_objects_list_store is created in the VismolGLWidget 
    #         file and does not depend on the maintreeview of the main window.
    #     """
    #     if vismol_object.type == "molecule":
    #         i = vismol_object.index 
    #         data = [vismol_object.active           , 
    #                 str(i)                  ,
    #                 vismol_object.name             , 
    #                 str(len(vismol_object.atoms))  , 
    #                 str(len(vismol_object.frames)) ,
    #                 ]
    #         self.vm_objects_list_store.append(data)
    
    def delete_by_index(self, index = None):
        """ Function doc """
        self.viewing_selections = []
        self.picking_selections = [None]*4
        self.vm_objects.pop(index)
        #self.vm_widget.updateGL()
    #""" 
    
    def select (self, vismol_object=None, indexes=None, sele=None):
        """ Function doc """
        if vismol_object is None:
            vismol_object = self.vm_objects[-1]
        
        if sele is None:
            sele = self.current_selection
        
        if indexes == "all":
            self.selections[sele].selecting_by_indexes(vismol_object=vismol_object,
                                                       indexes=range(0, int(len(vismol_object.atoms)/2)))
        self.vm_glcore.queue_draw()
        
    def orient (self, obj =  None):
        """ Function doc """  
    
    def center (self, vismol_object):
        """ Function doc """
        print ("center", vismol_object)
        frame = self.get_frame ()
        vismol_object._get_center_of_mass (frame)
        self.vm_glcore.center_on_coordinates(vismol_object, vismol_object.mass_center)


    def center_by_index(self, Vobject =  None, index = None):
        """ Function doc """  
        #mass_center = self.vm_objects[index].mass_center
        #self.vm_widget.center_on_atom(mass_center)
        mass_center = self.vm_objects_dic[index].mass_center
        
    def disable_by_index (self, index = 0 ):#, dictionary = False):
        """When the variable "dictionary" is active, the function accesses 
        a vismol object through the dictionary "self.vm_objects_dic". 
        Each vismol object has a unique access key (int), which, in 
        easyhybrid, is generated in the method: add_vismol_object_to_vismol_session.

        In the vismol interface the enable_by_index/disable_by_index methods
        access the vismol objects by their position in the "self.vm_objects" 
        list (this is because when an object is deleted in the vismol 
        interface, the treeview"s liststore is rewritten) """
        #if dictionary:
        #    self.vm_objects_dic[index].active = False
        #else:
        #    self.vm_objects[index].active = False
        #self.vm_widget.queue_draw()
        
        self.vm_objects_dic[index].active = False
        self.vm_glcore.queue_draw()

    def enable_by_index (self, index = 0):#, dictionary = True):
        """When the variable "dictionary" is active, the function accesses 
        a vismol object through the dictionary "self.vm_objects_dic". 
        Each vismol object has a unique access key (int), which, in 
        easyhybrid, is generated in the method: add_vismol_object_to_vismol_session.

        In the vismol interface the enable_by_index/disable_by_index methods
        access the vismol objects by their position in the "self.vm_objects" 
        list (this is because when an object is deleted in the vismol 
        interface, the treeview"s liststore is rewritten) """
        
        #if dictionary:
        #    self.vm_objects_dic[index].active = True
        #else:
        #    self.vm_objects[index].active = True
        #self.vm_widget.queue_draw()
        self.vm_objects_dic[index].active = True
        self.vm_glcore.queue_draw()
        
        
    def edit_by_index(self, index = 0):
        """ Function doc """
        #self.vm_objects[index].editing = not self.vm_objects[index].editing
        self.vm_objects_dic[index].editing = not self.vm_objects_dic[index].editing
        #self.vm_widget.queue_draw()
    
    def set_color_by_index (self, vismol_object = None, indexes = [ ], color = [0.9, 0.9, 0.9] ):
        """ Function doc """
        #selection         = self.selections[self.current_selection]
        
        #fixedlist = []
        #if len(indexes) > 0:
        
        for atom_index in indexes:
            vismol_object.atoms[atom_index].color = color    
            print(atom_index, color)
        print(vismol_object.colors)
        vismol_object._generate_color_vectors ( do_colors         = True,
                                                do_colors_idx     = False,
                                                do_colors_raindow = False,
                                                do_vdw_dot_sizes  = False,
                                                do_cov_dot_sizes  = False,
                                               )
        print(vismol_object.colors)

        self.vm_glcore.queue_draw()
        for rep  in vismol_object.representations.keys():
            if vismol_object.representations[rep]:
                try:
                    vismol_object.representations[rep]._set_colors_to_buffer()
                except:
                    print('"VisMol/vModel/Representations.py, line 123, in _set_colors_to_buffer GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ctypes.ArgumentError: argument 2: <class "TypeError">: wrong type"')
                    
        return True


        #refresh = self.main_session.pDynamo_session.define_free_or_fixed_atoms_from_iterable (fixedlist)
    
    def set_frame (self, frame = 0):
        """ Function doc """
        self.vm_glcore.frame = frame
        self.vm_glcore.queue_draw()

        #self.vm_widget.updateGL()
    
    def get_distance(self):
        """ Function doc """
        if self._picking_selection_mode:
            print(self.picking_selections.picking_selections_list)
    
    def get_frame (self):
        """ Function doc """
        #""" Function doc """
        frame = self.vm_glcore.frame
        return frame
    
    def get_vismol_object_dict (self):
        """ Function doc """
        Vobjects_dic = {}
    
        for vobj_id, Vobject in self.vm_objects_dic.items():
            #print ("----------------------- > get_vismol_object_list ", Vobject.label)
            index = self.vm_objects.index(Vobject)
            name = Vobject.label
            ##print( "\n label get_vismol_object_list:", name, index, len(Vobject.atoms) )
            Vobjects_dic[index] = name
    
        return Vobjects_dic
   
    def viewing_selection_mode(self, sel_type = "atom"):
        """ Function doc """        
        
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

       
