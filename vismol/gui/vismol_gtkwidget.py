#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_gtkwidget.py
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

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from libgl.vismol_glcore import VismolGLCore


class VismolGTKWidget(Gtk.GLArea):
    """ Object that contains the GLArea from GTK3+.
        It needs a vertex and shader to be created, maybe later I"ll
        add a function to change the shaders.
    """
    
    def __init__(self, vismol_session=None, width=640.0, height=420.0):
        """ Class initialiser
        """
        super(VismolGTKWidget, self).__init__()
        self.connect("realize", self.initialize)
        self.connect("render", self.render)
        self.connect("resize", self.reshape)
        self.connect("key-press-event", self.key_pressed)
        self.connect("key-release-event", self.key_released)
        self.connect("button-press-event", self.mouse_pressed)
        self.connect("button-release-event", self.mouse_released)
        self.connect("motion-notify-event", self.mouse_motion)
        self.connect("scroll-event", self.mouse_scroll)
        self.grab_focus()
        self.set_events(self.get_events() | Gdk.EventMask.SCROLL_MASK
                        | Gdk.EventMask.BUTTON_PRESS_MASK | Gdk.EventMask.BUTTON_RELEASE_MASK
                        | Gdk.EventMask.POINTER_MOTION_MASK | Gdk.EventMask.POINTER_MOTION_HINT_MASK
                        | Gdk.EventMask.KEY_PRESS_MASK | Gdk.EventMask.KEY_RELEASE_MASK)
        # self.vm_objects_list_store = Gtk.ListStore(bool,  # visible? 
        #                                            str,  # id
        #                                            str,  # name
        #                                            str,  # num of atoms
        #                                            str)  # num of frames
        self.vm_selection_modes_list_store = Gtk.ListStore(str)
        self.vm_session = vismol_session
        self.vm_glcore = VismolGLCore(self, vismol_session, width, height)
        self.glMenu_sele = None
        self.glMenu_bg = None
        self.glMenu_obj = None
    
    def initialize(self, widget):
        """ Enables the buffers and other charasteristics of the OpenGL context.
            sets the initial projection and view matrix
            
            self.flag -- Needed to only create one OpenGL program, otherwise a bunch of
                         programs will be created and use system resources. If the OpenGL
                         program will be changed change this value to True
        """
        if self.get_error() != None:
            print(self.get_error().args)
            print(self.get_error().code)
            print(self.get_error().domain)
            print(self.get_error().message)
            Gtk.main_quit()
        self.vm_glcore.initialize()
    
    def reshape(self, widget, width, height):
        """ Resizing function, takes the widht and height of the widget
            and modifies the view in the camera acording to the new values
        
            Keyword arguments:
            widget -- The widget that is performing resizing
            width -- Actual width of the window
            height -- Actual height of the window
        """
        self.vm_glcore.resize_window(width, height)
        self.queue_draw()
    
    def render(self, area, context):
        """ This is the function that will be called everytime the window
            needs to be re-drawed.
        """
        self.vm_glcore.render()
    
    def key_pressed(self, widget, event):
        """ The key_pressed function serves, as the names states, to catch
            events in the keyboard, e.g. letter "l" pressed, "backslash"
            pressed. Note that there is a difference between "A" and "a".
            Here I use a specific handler for each key pressed after
            discarding the CONTROL, ALT and SHIFT keys pressed (usefull
            for customized actions) and maintained, i.e. it"s the same as
            using Ctrl+Z to undo an action.
        """
        try:
            func = getattr(self, "_pressed_" + Gdk.keyval_name(event.keyval))
            func()
        except AttributeError as ae:
            pass
            # print(ae)
    
    def key_released(self, widget, event):
        """ Used to indicates a key has been released.
        """
        try:
            func = getattr(self, "_released_" + Gdk.keyval_name(event.keyval))
            func()
        except AttributeError as ae:
            pass
            # print(ae)
    
    def _pressed_Control_L(self):
        """ Function doc
        """
        self.vm_glcore.ctrl = True
    
    def _released_Control_L(self):
        """ Function doc
        """
        self.vm_glcore.ctrl = False
    
    def _pressed_Shift_L(self):
        """ Function doc
        """
        self.vm_glcore.shift = True
    
    def _released_Shift_L(self):
        """ Function doc
        """
        self.vm_glcore.shift = False
    
    def mouse_pressed(self, widget, event):
        """ Function doc
        """
        self.vm_glcore.mouse_pressed(event.button, event.x, event.y)
    
    def mouse_released(self, widget, event):
        """ Function doc
        """
        self.vm_glcore.mouse_released(event.button, event.x, event.y)
    
    def mouse_motion(self, widget, event):
        """ Function doc
        """
        self.vm_glcore.mouse_motion(event.x, event.y)
    
    def mouse_scroll(self, widget, event):
        """ Function doc
        """
        if event.direction == Gdk.ScrollDirection.UP:
            self.vm_glcore.mouse_scroll(1)
        if event.direction == Gdk.ScrollDirection.DOWN:
            self.vm_glcore.mouse_scroll(-1)
    
    
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
        self.build_glmenu(bg_menu=bg_menu, sele_menu=sele_menu,
                                    obj_menu=obj_menu, pick_menu=pick_menu)
    

    def build_submenus_from_dicts(self, menu_dict):
        """ Function doc
        """
        menu = Gtk.Menu()
        
        for key in menu_dict:
            mitem = Gtk.MenuItem(key)
            
            if menu_dict[key][0] == "submenu":
                #print(key)
                menu2 = self.build_submenus_from_dicts (menu_dict[key][1])
                mitem.set_submenu(menu2)
            
            
            elif menu_dict[key][0] == "separator":
                mitem = Gtk.SeparatorMenuItem()
                #menu2 = self.build_submenus_from_dicts (menu_dict[key][1])
                #mitem.set_submenu(menu2)
                #print(key)

            
            else:
                if menu_dict[key][1] != None:
                    mitem.connect("activate", menu_dict[key][1])
                else:
                    pass
            menu.append(mitem)
        
        return menu
        #menu.show_all()

    def build_glmenu_from_dicts (self, menu_dict, glMenu):
        """ Function doc """
        for key in menu_dict:
            mitem = Gtk.MenuItem(label = key)
            
            if menu_dict[key][0] == "submenu":
                menu2 = self.build_submenus_from_dicts (menu_dict[key][1])
                mitem.set_submenu(menu2)
            
            elif menu_dict[key][0] == "separator":
                mitem = Gtk.SeparatorMenuItem()
          
            else:
                if menu_dict[key][1] != None:
                    mitem.connect("activate", menu_dict[key][1])
                else:
                    pass
            glMenu.append(mitem) 

    def build_glmenu (self,  bg_menu  = None, sele_menu = None, obj_menu = None , pick_menu =  None):
        """ Function doc """
        
        """ Selection Menu """
        # --------------------------------------------------------------- #
        if sele_menu:
            self.glMenu_sele           = Gtk.Menu()
            self.glMenu_sele_toplabel =  Gtk.MenuItem(label = "selection")
            self.glMenu_sele.append (self.glMenu_sele_toplabel)
            
            self.build_glmenu_from_dicts( sele_menu, self.glMenu_sele)
           
            self.glMenu_sele.show_all()

        else:
            self.glMenu_sele = None
        # --------------------------------------------------------------- #
        
        """ Picking Menu """
        # --------------------------------------------------------------- #
        if pick_menu:
            self.glMenu_pick           = Gtk.Menu()
            self.glMenu_pick_toplabel =  Gtk.MenuItem(label = "picking")
            self.glMenu_pick.append (self.glMenu_pick_toplabel)
            
            self.build_glmenu_from_dicts( pick_menu, self.glMenu_pick)
           
            self.glMenu_pick.show_all()

        else:
            self.glMenu_pick = None
        # --------------------------------------------------------------- #

        """ Background Menu """
        # --------------------------------------------------------------- #
        if bg_menu:
            self.glMenu_bg  = Gtk.Menu()
            self.glMenu_bg_toplabel =  Gtk.MenuItem(label = "background")
            self.glMenu_bg.append (self.glMenu_bg_toplabel)

            self.build_glmenu_from_dicts( bg_menu, self.glMenu_bg)
            #self.glMenu_bg = self.build_submenus_from_dicts (bg_menu)
            self.glMenu_bg.show_all()
        else:
            self.glMenu_bg = None
        
        
        if obj_menu:
            self.glMenu_obj  = Gtk.Menu()
            self.glMenu_obj_toplabel =  Gtk.MenuItem(label = "atom")
            self.glMenu_obj.append (self.glMenu_obj_toplabel)

            self.build_glmenu_from_dicts( obj_menu, self.glMenu_obj)
            self.glMenu_obj.show_all()
        else:
            self.glMenu_obj = None
            
     
    def show_gl_menu (self, signals = None, menu_type = None, info = None):
        """ Function doc """
        
        if menu_type == "bg_menu":
            if self.glMenu_bg:
                self.glMenu_bg.popup(None, None, None, None, 0, 0)
        
        
        if menu_type == "sele_menu":
            if self.glMenu_sele:
                self.glMenu_sele.popup(None, None, None, None, 0, 0)
        
        if menu_type == "pick_menu":
            if self.glMenu_pick:
                self.glMenu_pick.popup(None, None, None, None, 0, 0)
        
        
        if menu_type == "obj_menu":
            if self.glMenu_obj:
                self.glMenu_obj_toplabel.set_label(info)
                self.glMenu_obj.popup(None, None, None, None, 0, 0)
        
    