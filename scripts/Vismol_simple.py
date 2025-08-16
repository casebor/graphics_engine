#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import gi
import sys
import logging
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from vismol.libgl.vismol_gl_main import VismolGLMain
from vismol.gui.vismol_widget_main import VismolWidgetMain
from vismol.gui.vismol_widget_builder import VismolWidgetBuilder
from vismol.core.vismol_config import VismolConfig
from vismol.core.vismol_session_main import VismolSessionMain
from vismol.core.vismol_session_builder import VismolSessionBuilder


logger = logging.getLogger(__name__)


class MainVismolWindow(Gtk.Window):
    """docstring for MainVismolWindow"""
    def __init__(self):
        super(MainVismolWindow, self).__init__(title="Vismol Main")
        self.set_default_size(640, 480)
        self.main_container = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.add(self.main_container)
        
        self.main_config = VismolConfig()
        logger.setLevel(self.main_config.console_log_level)
        self.main_widget = VismolWidgetMain(vismol_config=self.main_config, width=640, height=420)
        self.main_session = VismolSessionMain(vismol_widget=self.main_widget, vismol_config=self.main_config)
        self.main_widget.vm_session = self.main_session
        self.main_widget.vm_glcore.vm_session = self.main_session
        self.main_widget.insert_glmenu()
        self.main_container.pack_start(self.main_widget, True, True, 0)
        
        self.button_container = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        
        self.button = Gtk.Button(label="Open Builder Window")
        self.button.connect("clicked", self.on_button_clicked)
        self.button.show()
        self.button_container.pack_start(self.button, True, True, 0)
        self.button_container.show()
        self.main_container.pack_start(self.button_container, False, False, 0)
        self.main_container.show()
        self.connect("key-press-event", self.main_widget.key_pressed)
        self.connect("key-release-event", self.main_widget.key_released)
        self.connect("delete-event", Gtk.main_quit)
    
    def on_button_clicked(self, widget):
        builder_window = BuilderVismolWindow()
        builder_window.show()

class BuilderVismolWindow(Gtk.Window):
    """docstring for BuilderVismolWindow"""
    def __init__(self):
        super(BuilderVismolWindow, self).__init__(title="Builder Window")
        self.set_default_size(420, 360)
        my_container = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.add(my_container)
        
        builder_config = VismolConfig(builder=True)
        builder_widget = VismolWidgetBuilder(builder_config, width=420, height=360)
        builder_session = VismolSessionBuilder(vismol_widget=builder_widget, vismol_config=builder_config)
        builder_widget.vm_session = builder_session
        builder_widget.vm_glcore.vm_session = builder_session
        my_container.pack_start(builder_widget, True, True, 0)
        my_container.show()
        self.connect("key-press-event", builder_widget.key_pressed)
        self.connect("key-release-event", builder_widget.key_released)
        self.connect("destroy", self.on_destroy)
    
    def on_destroy(self, widget):
        pass

def main():
    # logging.basicConfig(format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
    #                     datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)
    logging.basicConfig(format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
                        datefmt="%H:%M:%S", level=logging.DEBUG)
    multiple_window = MainVismolWindow()
    multiple_window.show()
    try:
        filein = sys.argv[-1]
        multiple_window.vm_session.load_molecule(filein)
    except:
        pass
    Gtk.main()
    return 0

if __name__ == "__main__":
    main()
