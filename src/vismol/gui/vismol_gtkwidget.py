#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from logging import getLogger


logger = getLogger(__name__)


class VismolGTKWidget(Gtk.GLArea):
    """ Object that contains the GLArea from GTK3+.
        It needs a vertex and shader to be created, maybe later I"ll
        add a function to change the shaders.
    """
    
    def __init__(self, vismol_config, width, height):
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
        self.show()
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
        # self.vm_selection_modes_list_store = Gtk.ListStore(str)
        # self.vm_session = vismol_session
        # self.vm_glcore = VismolGLMain(self, vismol_session, width, height)
        # self.glmenu_bg = None
        # self.glmenu_sele = None
        # self.glmenu_obj = None
        # self.glmenu_pick = None
        # self.filechooser = None
        # self.selection_box_frame = None

    def initialize(self, widget):
        logger.critical("NotImplementedError, the child class must implement initialize")
        raise NotImplementedError("Subclasses must implement this method")
    
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
    
    def mouse_pressed(self, widget, event):
        """ Function doc """
        self.vm_glcore.mouse_pressed(event.button, event.x, event.y)
    
    def mouse_released(self, widget, event):
        """ Function doc """
        self.vm_glcore.mouse_released(event.button, event.x, event.y)
    
    def mouse_motion(self, widget, event):
        """ Function doc """
        self.vm_glcore.mouse_motion(event.x, event.y)
    
    def mouse_scroll(self, widget, event):
        """ Function doc
        """
        if event.direction == Gdk.ScrollDirection.UP:
            self.vm_glcore.mouse_scroll(1)
        if event.direction == Gdk.ScrollDirection.DOWN:
            self.vm_glcore.mouse_scroll(-1)
    
    def open_file(self, widget):
        logger.critical("NotImplementedError, the child class must implement open_file")
        raise NotImplementedError("Subclasses must implement this method")
    
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
            logger.debug("Press key {} has not been assigned to a handler "\
                         "yet".format(Gdk.keyval_name(event.keyval)))
            logger.error(ae)

    def key_released(self, widget, event):
        """ Used to indicates a key has been released.
        """
        try:
            func = getattr(self, "_released_" + Gdk.keyval_name(event.keyval))
            func()
        except AttributeError as ae:
            logger.debug("Release key {} has not been assigned to a handler "\
                         "yet".format(Gdk.keyval_name(event.keyval)))
            logger.error(ae)
    
    def _pressed_Escape(self):
        logger.critical("NotImplementedError, the child class must implement _pressed_Escape")
        raise NotImplementedError("Subclasses must implement this method")
    
    def quit(self):
        logger.critical("NotImplementedError, the child class must implement quit")
        raise NotImplementedError("Subclasses must implement this method")
    
