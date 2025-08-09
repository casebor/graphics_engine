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
        super(MainVismolWindow, self).__init__(title="Vismol Window")
        self.set_default_size(640, 480)
        self.main_container = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.add(self.main_container)
        
        self.main_config = VismolConfig()
        self.main_widget = VismolWidgetMain(vismol_config=self.main_config, width=640.0, height=420.0)
        # self.main_widget.set_size_request(640, 400)
        self.main_session = VismolSessionMain(vismol_widget=self.main_widget, vismol_config=self.main_config)
        self.main_widget.vm_session = self.main_session
        self.main_widget.vm_glcore.vm_session = self.main_session
        self.main_widget.insert_glmenu()
        self.main_container.pack_start(self.main_widget, True, True, 0)
        
        self.button_container = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        
        self.button = Gtk.Button(label="Open Secondary Window")
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
        edit_window = SecondaryVismolWindow()
        edit_window.show()

class SecondaryVismolWindow(Gtk.Window):
    """docstring for SecondaryVismolWindow"""
    def __init__(self):
        super(SecondaryVismolWindow, self).__init__(title="Builder Window")
        self.set_default_size(420, 360)
        my_container = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.add(my_container)
        
        builder_config = VismolConfig(builder=True)
        builder_widget = VismolWidgetBuilder(builder_config, width=420.0, height=360.0)
        builder_session = VismolSessionBuilder(vismol_widget=builder_widget, vismol_config=builder_config)
        builder_widget.vm_session = builder_session
        builder_widget.vm_glcore.vm_session = builder_session
        my_container.pack_start(builder_widget, True, True, 0)
        my_container.show()
        self.connect("key-press-event", builder_widget.key_pressed)
        self.connect("key-release-event", builder_widget.key_released)
        self.connect("destroy", self.on_destroy)
    
    def on_button_clicked(self, widget):
        print("Button pressed in secondary window")
    
    def on_destroy(self, widget):
        pass

def main():
    logging.basicConfig(format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
                        datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)
    # vm_session = VismolSession(toolkit="Gtk_3.0")
    # vm_session.vm_widget.insert_glmenu()
    # main_window = Gtk.Window(title="Vismol Window")
    # main_container = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    # main_container.pack_start(vm_session.vm_widget, True, True, 0)
    # main_container.show()
    # test_pop = TestPopUp()
    # main_container.pack_start(test_pop.vbox, True, True, 0)
    # build_container = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    # build_container.pack_start(test_pop.vbox, True, True, 0)
    # build_container.show()
    # main_window.connect("key-press-event", vm_session.vm_widget.key_pressed)
    # main_window.connect("key-release-event", vm_session.vm_widget.key_released)
    # main_window.add(main_container)
    # # main_window.add(test_pop.vbox)
    # main_window.connect("delete-event", Gtk.main_quit)
    # main_window.show()
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
