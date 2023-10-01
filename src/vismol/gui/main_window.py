#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  

import os
import gi, sys
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
import vismol.gui.icons as vg_ico
import vismol.gui.gtk_widgets as vg_gtk


class VismolMainWindow():
    """ Class doc """
    def __init__ (self, vismol_session=None):
        """ Class initialiser """
        self.vm_session = vismol_session
        self.vm_session.vm_widget.insert_glmenu()
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(os.path.dirname(vg_gtk.__file__), "MainWindow.glade"))
        self.builder.connect_signals(self)
        self.window = self.builder.get_object("MainWindow")
        self.window.set_title("Vismol 1.0")
        self.window.set_default_size(800, 600)
        # self.gtkbox = self.builder.get_object("MainBox")
        # self.gtkbox.pack_start(self.vm_session.vm_widget, True, True, 0)
        iconpath = os.path.dirname(vg_ico.__file__)
        self.main_iconbar = self.builder.get_object("MainIconBar")
        
        self.main_box = self.builder.get_object("MainBox")
        self.main_box.pack_start(self.vm_session.vm_widget, True, True, 0)
        # self.main_statusbar = self.builder.get_object("MainStatusBar")
        # self.main_statusbar.push(1, "Welcome to Vismol")
        
        self.window.connect("key-press-event", self.vm_session.vm_widget.key_pressed)
        self.window.connect("key-release-event", self.vm_session.vm_widget.key_released)
        self.window.connect("delete-event", Gtk.main_quit)
        self.window.show_all()