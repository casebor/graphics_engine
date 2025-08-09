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
    
    def __init__(self, vismol_config, width, height):
        """ Class initialiser
        """
        super(VismolWidgetBuilder, self).__init__(vismol_config, width, height)
        self.vm_glcore = VismolGLBuilder(self, vismol_config, width, height)
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
        self.vm_session = None
    
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
        self.vm_glcore.initialize_builder()
    
    def _pressed_p(self):
        """ Function doc """
        for atom in self.vm_session.vm_objects_dic[0].molecule.atoms.values():
            print(atom, atom.coords())
    