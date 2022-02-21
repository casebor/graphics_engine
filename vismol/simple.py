#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  simple.py
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

import gi, sys
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
from core.vismol_session import VismolSession


def main():
    vm_session = VismolSession(glwidget=True, toolkit="gtk3")
    vm_session.insert_glmenu()
    window = Gtk.Window(title="Vismol window")
    container = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    container.pack_start(vm_session.vm_widget, True, True, 0)
    window.connect("key-press-event", vm_session.vm_widget.key_pressed)
    window.connect("key-release-event", vm_session.vm_widget.key_released)
    window.add(container)
    window.connect("delete-event", Gtk.main_quit)
    window.show_all()
    try:
        filein = sys.argv[-1]
        vm_session.load(filein)
    except:
        pass
    Gtk.main()
    return 0

if __name__ == '__main__':
    main()
