#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from logging import getLogger
from vismol.gui.filechooser import FileChooser
from vismol.libgl.vismol_gl_main import VismolGLMain
from vismol.gui.vismol_gtkwidget import VismolGTKWidget


logger = getLogger(__name__)


class VismolWidgetMain(VismolGTKWidget):
    """ Object that contains the GLArea from GTK3+.
        It needs a vertex and shader to be created, maybe later I"ll
        add a function to change the shaders.
    """
    
    def __init__(self, vismol_config, width=640.0, height=420.0):
        """ Class initialiser
        """
        super(VismolWidgetMain, self).__init__(vismol_config, width, height)
        self.vm_glcore = VismolGLMain(self, vismol_config, width, height)
        # self.connect("realize", self.initialize)
        # self.connect("render", self.render)
        # self.connect("resize", self.reshape)
        # self.connect("key-press-event", self.key_pressed)
        # self.connect("key-release-event", self.key_released)
        # self.connect("button-press-event", self.mouse_pressed)
        # self.connect("button-release-event", self.mouse_released)
        # self.connect("motion-notify-event", self.mouse_motion)
        # self.connect("scroll-event", self.mouse_scroll)
        # self.show()
        # self.grab_focus()
        # self.set_events(self.get_events() | Gdk.EventMask.SCROLL_MASK
        #                 | Gdk.EventMask.BUTTON_PRESS_MASK | Gdk.EventMask.BUTTON_RELEASE_MASK
        #                 | Gdk.EventMask.POINTER_MOTION_MASK | Gdk.EventMask.POINTER_MOTION_HINT_MASK
        #                 | Gdk.EventMask.KEY_PRESS_MASK | Gdk.EventMask.KEY_RELEASE_MASK)
        # self.vm_objects_list_store = Gtk.ListStore(bool,  # visible? 
        #                                            str,  # id
        #                                            str,  # name
        #                                            str,  # num of atoms
        #                                            str)  # num of frames
        self.vm_selection_modes_list_store = Gtk.ListStore(str)
        self.vm_session = None
        # self.vm_glcore = VismolGLMain(self, vismol_session, width, height)
        self.glmenu_bg = None
        self.glmenu_sele = None
        self.glmenu_obj = None
        self.glmenu_pick = None
        self.filechooser = None
        self.selection_box_frame = None
    
    def initialize(self, widget):
        """ Enables the buffers and other charasteristics of the OpenGL context.
            sets the initial projection and view matrix
            
        """
        if self.get_error() != None:
            logger.critical(self.get_error().args)
            logger.critical(self.get_error().code)
            logger.critical(self.get_error().domain)
            logger.critical(self.get_error().message)
            Gtk.main_quit()
        self.vm_glcore.initialize_main()
    
    # def reshape(self, widget, width, height):
    #     """ Resizing function, takes the widht and height of the widget
    #         and modifies the view in the camera acording to the new values
        
    #         Keyword arguments:
    #         widget -- The widget that is performing resizing
    #         width -- Actual width of the window
    #         height -- Actual height of the window
    #     """
    #     self.vm_glcore.resize_window(width, height)
    #     self.queue_draw()
    
    # def render(self, area, context):
    #     """ This is the function that will be called everytime the window
    #         needs to be re-drawed.
    #     """
    #     self.vm_glcore.render()
    
    # def mouse_pressed(self, widget, event):
    #     """ Function doc """
    #     self.vm_glcore.mouse_pressed(event.button, event.x, event.y)
    
    # def mouse_released(self, widget, event):
    #     """ Function doc """
    #     self.vm_glcore.mouse_released(event.button, event.x, event.y)
    
    # def mouse_motion(self, widget, event):
    #     """ Function doc """
    #     self.vm_glcore.mouse_motion(event.x, event.y)
    
    # def mouse_scroll(self, widget, event):
    #     """ Function doc
    #     """
    #     if event.direction == Gdk.ScrollDirection.UP:
    #         self.vm_glcore.mouse_scroll(1)
    #     if event.direction == Gdk.ScrollDirection.DOWN:
    #         self.vm_glcore.mouse_scroll(-1)
    
    def _build_glmenu(self, bg_menu=None, sele_menu=None, obj_menu=None, pick_menu=None):
        """ Function doc """
        if bg_menu is not None:
            self.glmenu_bg = Gtk.Menu()
            self.glmenu_bg_toplabel = Gtk.MenuItem(label="background")
            self._build_glmenu_from_dicts(bg_menu, self.glmenu_bg)
            self.glmenu_bg.show_all()
        else:
            self.glmenu_bg = None
        
        if sele_menu is not None:
            self.glmenu_sele = Gtk.Menu()
            self.glmenu_sele_toplabel = Gtk.MenuItem(label="selection")
            self._build_glmenu_from_dicts(sele_menu, self.glmenu_sele)
            self.glmenu_sele.show_all()
        else:
            self.glmenu_sele = None
        
        if pick_menu is not None:
            self.glmenu_pick = Gtk.Menu()
            self.glmenu_pick_toplabel = Gtk.MenuItem(label="picking")
            self.glmenu_pick.append(self.glmenu_pick_toplabel)
            self._build_glmenu_from_dicts(pick_menu, self.glmenu_pick)
            self.glmenu_pick.show_all()
        else:
            self.glmenu_pick = None
        
        if obj_menu is not None:
            self.glmenu_obj = Gtk.Menu()
            self.glmenu_obj_toplabel = Gtk.MenuItem(label="object")
            self.glmenu_obj.append(self.glmenu_obj_toplabel)
            self._build_glmenu_from_dicts(obj_menu, self.glmenu_obj)
            self.glmenu_obj.show_all()
        else:
            self.glmenu_obj = None
     
    def open_file(self, widget):
        """ Function doc """
        if self.filechooser is None:
            self.filechooser = FileChooser()
        filename = self.filechooser.open()
        self.vm_session.load_molecule(filename)
    
    # def key_pressed(self, widget, event):
    #     """ The key_pressed function serves, as the names states, to catch
    #         events in the keyboard, e.g. letter "l" pressed, "backslash"
    #         pressed. Note that there is a difference between "A" and "a".
    #         Here I use a specific handler for each key pressed after
    #         discarding the CONTROL, ALT and SHIFT keys pressed (usefull
    #         for customized actions) and maintained, i.e. it"s the same as
    #         using Ctrl+Z to undo an action.
    #     """
    #     try:
    #         func = getattr(self, "_pressed_" + Gdk.keyval_name(event.keyval))
    #         func()
    #     except AttributeError as ae:
    #         logger.debug("Press key {} has not been assigned to a handler "\
    #                      "yet".format(Gdk.keyval_name(event.keyval)))
    
    # def key_released(self, widget, event):
    #     """ Used to indicates a key has been released.
    #     """
    #     try:
    #         func = getattr(self, "_released_" + Gdk.keyval_name(event.keyval))
    #         func()
    #     except AttributeError as ae:
    #         logger.debug("Release key {} has not been assigned to a handler "\
    #                      "yet".format(Gdk.keyval_name(event.keyval)))
    
    def _pressed_p(self):
        """ Function doc """
        for vm_object in self.vm_session.vm_objects_dic.values():
            for i, atom in vm_object.molecule.atoms.items():
                print(i, atom)
    
    def _pressed_Escape(self):
        """ Function doc """
        self.quit()
    
    def _pressed_Right(self):
        """ Function doc """
        self.vm_session.forward_frame()
        self.queue_draw()
    
    def _released_Right(self):
        """ Function doc """
        pass
    
    def _pressed_Left(self):
        """ Function doc """
        self.vm_session.reverse_frame()
        self.queue_draw()
    
    def _released_Left(self):
        """ Function doc """
        pass
    
    def _pressed_Control_L(self):
        """ Function doc """
        self.vm_glcore.ctrl = True
    
    def _released_Control_L(self):
        """ Function doc """
        self.vm_glcore.ctrl = False
    
    def _pressed_Shift_L(self):
        """ Function doc """
        self.vm_glcore.shift = True
    
    def _released_Shift_L(self):
        """ Function doc """
        self.vm_glcore.shift = False
    
    def _selection_type_picking(self, widget):
        if self.selection_box_frame:
            self.selection_box_frame.change_toggle_button_selecting_mode_status(True)
        else:
            self.vm_session.picking_selection_mode = True
        self.queue_draw()
    
    def _selection_type_viewing(self, widget):
        if self.selection_box_frame:
            self.selection_box_frame.change_toggle_button_selecting_mode_status(False)
        else:
            self.vm_session.picking_selection_mode = False
        self.queue_draw()
    
    def quit(self):
        logger.info("Thank you for using our software :)")
        logger.info("Quitting Vismol")
        Gtk.main_quit()
    
    def _viewing_selection_mode_atom(self, widget):
        self.vm_session.viewing_selection_mode(sel_type="atom")
    
    def _viewing_selection_mode_residue(self, widget):
        self.vm_session.viewing_selection_mode(sel_type="residue")
    
    def _viewing_selection_mode_chain(self, widget):
        self.vm_session.viewing_selection_mode(sel_type="chain")
    
    def menu_show_dots(self, widget):
        self.vm_session.show_or_hide(rep_type="dots", show=True)
    
    def menu_hide_dots(self, widget):
        self.vm_session.show_or_hide(rep_type="dots", show=False)
    
    def menu_show_lines(self, widget):
        self.vm_session.show_or_hide(rep_type="lines", show=True)
    
    def menu_hide_lines(self, widget):
        self.vm_session.show_or_hide(rep_type="lines", show=False)
    
    def menu_show_nonbonded(self, widget):
        self.vm_session.show_or_hide(rep_type="nonbonded", show=True)
    
    def menu_hide_nonbonded(self, widget):
        self.vm_session.show_or_hide(rep_type="nonbonded", show=False)
    
    def menu_show_impostor(self, widget):
        self.vm_session.show_or_hide(rep_type="impostor", show=True)
    
    def menu_hide_impostor(self, widget):
        self.vm_session.show_or_hide(rep_type="impostor", show=False)
    
    def menu_show_spheres(self, widget):
        self.vm_session.show_or_hide(rep_type="spheres", show=True)
    
    def menu_hide_spheres(self, widget):
        self.vm_session.show_or_hide(rep_type="spheres", show=False)
    
    def menu_show_sticks(self, widget):
        self.vm_session.show_or_hide(rep_type="sticks", show=True)
    
    def menu_hide_sticks(self, widget):
        self.vm_session.show_or_hide(rep_type="sticks", show=False)
    
    def invert_selection(self, widget):
        self.vm_session.selections[self.vm_session.current_selection].invert_selection()
    
    def insert_glmenu(self, bg_menu=None, sele_menu=None, obj_menu=None, pick_menu=None):
        """ Function doc """
        if bg_menu is None:
            """ Standard Bg Menu"""
            bg_menu = {"Open File": ["MenuItem", self.open_file],
                       "separator": ["separator", None],
                       "Selection Mode": ["submenu",
                                            {"by atom": ["MenuItem", self._viewing_selection_mode_atom],
                                             "by residue": ["MenuItem", self._viewing_selection_mode_residue],
                                             "by chain": ["MenuItem", self._viewing_selection_mode_chain],
                                            }
                                         ],
                       "Selection Type": ["submenu",
                                            {"viewing": ["MenuItem", self._selection_type_viewing],
                                             "picking": ["MenuItem", self._selection_type_picking],
                                            }
                                         ],
                       "separator": ["separator", None],
                       "Quit": ["MenuItem", self.quit],
                      }
        
        if sele_menu is None:
            """ Standard Sele Menu """
            sele_menu = {"Show": ["submenu",
                                    {"dots": ["MenuItem", self.menu_show_dots],
                                     "lines": ["MenuItem", self.menu_show_lines],
                                     "nonbonded": ["MenuItem", self.menu_show_nonbonded],
                                     "sticks": ["MenuItem", self.menu_show_sticks],
                                     "impostor": ["MenuItem", self.menu_show_impostor],
                                     "spheres": ["MenuItem", self.menu_show_spheres],
                                    }
                                 ],
                         "Hide": ["submenu",
                                    {"dots": ["MenuItem", self.menu_hide_dots],
                                     "lines": ["MenuItem", self.menu_hide_lines],
                                     "nonbonded": ["MenuItem", self.menu_hide_nonbonded],
                                     "sticks": ["MenuItem", self.menu_hide_sticks],
                                     "impostor": ["MenuItem", self.menu_hide_impostor],
                                     "spheres": ["MenuItem", self.menu_hide_spheres],
                                    }
                                 ],
                         "separator":["separator", None],
                         "Invert Selection": ["MenuItem", self.invert_selection],
                        }
        
        if obj_menu is None:
            """ Standard Obj Menu"""
            obj_menu = {"Show": ["submenu",
                                    {"dots": ["MenuItem", None],
                                     "lines": ["MenuItem", None],
                                     "nonbonded": ["MenuItem", None],
                                    }
                                ],
                        "Hide": ["submenu",
                                    {"dots": ["MenuItem", None],
                                     "lines": ["MenuItem", None],
                                     "nonbonded": ["MenuItem", None],
                                    }
                                ],
                        "separator":["separator", None],
                        "label": ["submenu",
                                    {"Atom": ["submenu",
                                                {"index": ["MenuItem", None],
                                                 "name": ["MenuItem", None],
                                                 "residue": ["MenuItem", None],
                                                 "chain": ["MenuItem", None],
                                                 }
                                             ],
                                     "Residue": ["submenu",
                                                    {"index": ["MenuItem", None],
                                                     "name": ["MenuItem", None],
                                                     "chain": ["MenuItem", None],
                                                     }
                                                ],
                                     "Chain": ["submenu",
                                                {"name": ["MenuItem", None],
                                                 }
                                              ],
                                    },
                                 ],
                        }
        
        if pick_menu is None:
            """ Standard Sele Menu """
            pick_menu = {"Show": ["submenu",
                                    {"dots": ["MenuItem", None],
                                     "lines": ["MenuItem", None],
                                     "nonbonded": ["MenuItem", None],
                                    }
                                  ],
                         "Hide": ["submenu",
                                    {"dots": ["MenuItem", None],
                                     "lines": ["MenuItem", None],
                                     "nonbonded": ["MenuItem", None],
                                    }
                                  ],
                        }
        
        self._build_glmenu(bg_menu=bg_menu, sele_menu=sele_menu, obj_menu=obj_menu, pick_menu=pick_menu)
    
    def show_gl_menu(self, signals=None, menu_type=None, info=None):
        """ Function doc """
        if menu_type == "bg_menu":
            if self.glmenu_bg:
                self.glmenu_bg.popup(None, None, None, None, 0, 0)
        
        if menu_type == "sele_menu":
            if self.glmenu_sele:
                self.glmenu_sele.popup(None, None, None, None, 0, 0)
        
        if menu_type == "pick_menu":
            if self.glmenu_pick:
                self.glmenu_pick.popup(None, None, None, None, 0, 0)
        
        if menu_type == "obj_menu":
            if self.glmenu_obj:
                self.glmenu_obj_toplabel.set_label(info)
                self.glmenu_obj.popup(None, None, None, None, 0, 0)
    
    def _build_submenus_from_dicts(self, menu_dict):
        """ Function doc """
        menu = Gtk.Menu()
        for key in menu_dict:
            mitem = Gtk.MenuItem(key)
            if menu_dict[key][0] == "submenu":
                menu2 = self._build_submenus_from_dicts(menu_dict[key][1])
                mitem.set_submenu(menu2)
            elif menu_dict[key][0] == "separator":
                mitem = Gtk.SeparatorMenuItem()
            else:
                if menu_dict[key][1] != None:
                    mitem.connect("activate", menu_dict[key][1])
                else:
                    pass
            menu.append(mitem)
        return menu
    
    def _build_glmenu_from_dicts(self, menu_dict, glMenu):
        """ Function doc """
        for key in menu_dict:
            mitem = Gtk.MenuItem(label=key)
            if menu_dict[key][0] == "submenu":
                menu2 = self._build_submenus_from_dicts(menu_dict[key][1])
                mitem.set_submenu(menu2)
            elif menu_dict[key][0] == "separator":
                mitem = Gtk.SeparatorMenuItem()
            else:
                if menu_dict[key][1] != None:
                    mitem.connect("activate", menu_dict[key][1])
                else:
                    pass
            glMenu.append(mitem)
    
    def _pressed_Escape(self):
        """ Function doc """
        self.quit()
    
    def quit(self):
        logger.info("Thank you for using our software :)")
        logger.info("Quitting Vismol")
        Gtk.main_quit()
    
    def _released_p(self):
        """ Function doc """
        self.vm_session.vm_objects_dic[0].representations["spheres"].print_vbo_vao()
    
    def _released_c(self):
        """ Function doc """
        self.vm_glcore.glcamera.print_params()
        pass
    
