#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  __main__.py
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

import os
import gi, sys
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from vismol.gui.main_window import VismolMainWindow
from vismol.core.vismol_session import VismolSession


def main():
    """ Function doc """
    vm_session = VismolSession(toolkit="Gtk_3.0")
    vm_main = VismolMainWindow(vismol_session=vm_session)
    Gtk.main()
    return 0

if __name__ == "__main__":
    main()
