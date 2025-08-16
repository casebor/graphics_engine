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
    
    def __init__(self, vismol_config: "VismolConfig", width: int, height: int):
        """ Class initialiser
        """
        super(VismolGTKWidget, self).__init__()
        self.vm_config = vismol_config
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
    
    def initialize(self, widget: Gtk.GLArea) -> None:
        logger.critical("NotImplementedError, the child class must implement initialize")
        raise NotImplementedError("Subclasses must implement this method")
    
    def key_pressed(self, widget: Gtk.GLArea, event: Gdk.EventKey) -> None:
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
            logger.debug(ae)
    
    def key_released(self, widget: Gtk.GLArea, event: Gdk.EventKey) -> None:
        """ Used to indicates a key has been released.
        """
        try:
            func = getattr(self, "_released_" + Gdk.keyval_name(event.keyval))
            func()
        except AttributeError as ae:
            logger.debug("Release key {} has not been assigned to a handler "\
                         "yet".format(Gdk.keyval_name(event.keyval)))
            logger.debug(ae)
    
    def mouse_motion(self, widget: Gtk.GLArea, event: Gdk.EventMotion) -> None:
        """ Function doc """
        self.vm_glcore.mouse_motion(event.x, event.y)
    
    def mouse_pressed(self, widget: Gtk.GLArea, event: Gdk.EventButton) -> None:
        """ Function doc """
        self.vm_glcore.mouse_pressed(event.button, event.x, event.y)
    
    def mouse_released(self, widget: Gtk.GLArea, event: Gdk.EventButton) -> None:
        """ Function doc """
        self.vm_glcore.mouse_released(event.button, event.x, event.y)
    
    def mouse_scroll(self, widget: Gtk.GLArea, event: Gdk.EventScroll) -> None:
        """ Function doc
        """
        if event.direction == Gdk.ScrollDirection.UP:
            self.vm_glcore.mouse_scroll(1)
        if event.direction == Gdk.ScrollDirection.DOWN:
            self.vm_glcore.mouse_scroll(-1)
    
    def open_file(self, widget: Gtk.GLArea) -> None:
        logger.critical("NotImplementedError, the child class must implement open_file")
        raise NotImplementedError("Subclasses must implement this method")
    
    def render(self, area: Gtk.GLArea, context: gi.repository.GdkX11.X11GLContext) -> bool:
        """ This is the function that will be called everytime the window
            needs to be re-drawed.
        """
        self.vm_glcore.render()
    
    def reshape(self, widget: Gtk.GLArea, width: int, height: int) -> None:
        """ Resizing function, takes the widht and height of the widget
            and modifies the view in the camera acording to the new values
        
            Keyword arguments:
            widget -- The widget that is performing resizing
            width -- Actual width of the window
            height -- Actual height of the window
        """
        self.vm_glcore.resize_window(width, height)
        self.queue_draw()
    