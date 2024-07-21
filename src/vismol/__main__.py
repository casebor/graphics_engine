#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import gi, sys
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from vismol.gui.vismol_main import VismolMainWindow
from vismol.core.vismol_session import VismolSession


def main() -> int:
    """
    Main function for Vismol, has no arguments for now and uses gtk3 as default
    window manager.
    
    Returns:
        0: if there was no error.
    
    """
    vm_session = VismolSession(widget=True, toolkit="gtk3")
    if len(sys.argv) >= 2:
            filein = sys.argv[-1]
    else:
        filein = None
    vm_main = VismolMainWindow(vismol_session=vm_session, filein=filein)
    return 0

if __name__ == "__main__":
    main()
