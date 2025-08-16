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
    
    def __init__(self, vismol_config: "VismolConfig", width: int=640, height: int=420):
        """ Class initialiser
        """
        super(VismolWidgetMain, self).__init__(vismol_config, width, height)
        self.vm_glcore = VismolGLMain(self, vismol_config, width, height)
        self.vm_selection_modes_list_store = Gtk.ListStore(str)
        self.vm_session = None
        self.glmenu_bg = None
        self.glmenu_sele = None
        self.glmenu_obj = None
        self.glmenu_pick = None
        self.filechooser = None
        self.selection_box_frame = None
    
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
    
    def _build_glmenu_from_dicts(self, menu_dict, glMenu) -> None:
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
    
    def _pressed_Control_L(self) -> None:
        """ Function doc """
        self.vm_glcore.ctrl = True
    
    def _pressed_Escape(self) -> None:
        """ Function doc """
        self.quit()
    
    def _pressed_Left(self) -> None:
        """ Function doc """
        self.vm_session.reverse_frame()
        self.queue_draw()
    
    def _pressed_p(self) -> None:
        """ Function doc """
        # for vm_object in self.vm_session.vm_objects_dic.values():
        #     for i, atom in vm_object.molecule.atoms.items():
        #         print(i, atom)
        pass
    
    def _pressed_Right(self) -> None:
        """ Function doc """
        self.vm_session.forward_frame()
        self.queue_draw()
    
    def _pressed_Shift_L(self) -> None:
        """ Function doc """
        self.vm_glcore.shift = True
    
    def _released_c(self) -> None:
        """ Function doc """
        self.vm_glcore.glcamera.print_params()
    
    def _released_Control_L(self) -> None:
        """ Function doc """
        self.vm_glcore.ctrl = False
    
    def _released_Left(self) -> None:
        """ Function doc """
        pass
    
    def _released_p(self) -> None:
        """ Function doc """
        # self.vm_session.vm_objects_dic[0].representations["spheres"].print_vbo_vao()
        # for atom in self.vm_session.vm_objects_dic[0].molecule.atoms.values():
        #     print(atom)
        #     print("-"*60)
        #     for bond in atom.bonds:
        #         print(bond.atom_i, bond.atom_j)
        pass
    
    def _released_Right(self) -> None:
        """ Function doc """
        pass
    
    def _released_Shift_L(self) -> None:
        """ Function doc """
        self.vm_glcore.shift = False
    
    def _selection_type_picking(self, widget) -> None:
        if self.selection_box_frame:
            self.selection_box_frame.change_toggle_button_selecting_mode_status(True)
        else:
            self.vm_session.picking_selection_mode = True
        self.queue_draw()
    
    def _selection_type_viewing(self, widget) -> None:
        if self.selection_box_frame:
            self.selection_box_frame.change_toggle_button_selecting_mode_status(False)
        else:
            self.vm_session.picking_selection_mode = False
        self.queue_draw()
    
    def _viewing_selection_mode_atom(self, widget) -> None:
        self.vm_session.viewing_selection_mode(sel_type="atom")
    
    def _viewing_selection_mode_chain(self, widget) -> None:
        self.vm_session.viewing_selection_mode(sel_type="chain")
    
    def _viewing_selection_mode_residue(self, widget) -> None:
        self.vm_session.viewing_selection_mode(sel_type="residue")
    
    def initialize(self, widget: Gtk.GLArea) -> None:
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
    
    def insert_glmenu(self, bg_menu=None, sele_menu=None, obj_menu=None, pick_menu=None) -> None:
        """ Function doc """
        if bg_menu is None:
            """ Standard Bg Menu"""
            bg_menu = {"Open File": ["MenuItem", self.open_file],
                       "separator": ["separator", None],
                       "Selection Mode": [
                            "submenu",
                            {"by atom": ["MenuItem", self._viewing_selection_mode_atom],
                             "by residue": ["MenuItem", self._viewing_selection_mode_residue],
                             "by chain": ["MenuItem", self._viewing_selection_mode_chain],
                            }
                        ],
                       "Selection Type": [
                            "submenu",
                            {"viewing": ["MenuItem", self._selection_type_viewing],
                             "picking": ["MenuItem", self._selection_type_picking],
                            }
                        ],
                       "separator": ["separator", None],
                       "Quit": ["MenuItem", self.quit],
                      }
        
        if sele_menu is None:
            """ Standard Sele Menu """
            sele_menu = {"Show": [
                            "submenu",
                                {"dots": ["MenuItem", self.menu_show_dots],
                                 "lines": ["MenuItem", self.menu_show_lines],
                                 "nonbonded": ["MenuItem", self.menu_show_nonbonded],
                                 "sticks": ["MenuItem", self.menu_show_sticks],
                                 "impostor": ["MenuItem", self.menu_show_impostor],
                                 "spheres": ["MenuItem", self.menu_show_spheres],
                                }
                            ],
                         "Hide": [
                            "submenu",
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
    
    def invert_selection(self, widget) -> None:
        self.vm_session.selections[self.vm_session.current_selection].invert_selection()
    
    def open_file(self, widget) -> None:
        """ Function doc """
        if self.filechooser is None:
            self.filechooser = FileChooser()
        filename = self.filechooser.open()
        self.vm_session.load_molecule(filename)
    
    def menu_hide_dots(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="dots", show=False)
    
    def menu_hide_impostor(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="impostor", show=False)
    
    def menu_hide_lines(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="lines", show=False)
    
    def menu_hide_nonbonded(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="nonbonded", show=False)
    
    def menu_hide_spheres(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="spheres", show=False)
    
    def menu_hide_sticks(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="sticks", show=False)
    
    def menu_show_dots(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="dots", show=True)
    
    def menu_show_impostor(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="impostor", show=True)
    
    def menu_show_lines(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="lines", show=True)
    
    def menu_show_nonbonded(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="nonbonded", show=True)
    
    def menu_show_spheres(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="spheres", show=True)
    
    def menu_show_sticks(self, widget) -> None:
        self.vm_session.show_or_hide(rep_type="sticks", show=True)
    
    def show_gl_menu(self, signals=None, menu_type=None, info=None) -> None:
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
    
    def quit(self) -> None:
        logger.info("Thank you for using our software :)")
        logger.info("Quitting Vismol")
        Gtk.main_quit()
    