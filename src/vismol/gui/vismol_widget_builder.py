#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from logging import getLogger
from vismol.gui.vismol_gtkwidget import VismolGTKWidget
from vismol.libgl.vismol_gl_builder import VismolGLBuilder


logger = getLogger(__name__)


class VismolWidgetBuilder(VismolGTKWidget):
    """ Object that contains the GLArea from GTK3+.
        It needs a vertex and shader to be created, maybe later I"ll
        add a function to change the shaders.
    """
    
    def __init__(self, vismol_config: "VismolConfig", width: int, height: int):
        """ Class initialiser
        """
        super(VismolWidgetBuilder, self).__init__(vismol_config, width, height)
        self.vm_glcore = VismolGLBuilder(self, vismol_config, width, height)
        self.vm_session = None
    
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
        self.vm_glcore.initialize_builder()
    
    def _pressed_Escape(self, widget: Gtk.GLArea) -> None:
        """ Function doc """
        pass
    
    def _pressed_p(self) -> None:
        """ Function doc """
        for atom in self.vm_session.vm_objects_dic[0].molecule.atoms.values():
            print(atom, atom.coords())
    