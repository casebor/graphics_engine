#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

from sys import argv
import numpy as np
from libs import gtk_stuff
from libs.objects import GLContainer
import shaders.shaders as shaders

if __name__ == "__main__":
    """ Function doc """
    # if argv[1].endswith("pdb"):
    #     data = parse_pdb(argv[1])
    # elif argv[1].endswith("obj"):
    #     data = parse_obj(argv[1])
    # test = gtk_stuff.VMWindow(draw_type="points")
    test = gtk_stuff.VMWindow(vertex_shader=shaders.vertex_shader,
           fragment_shader = shaders.fragment_shader)
    # test = gtk_stuff.VMWindow(vertex_shader=shaders.vertex_shader_triangles,
    #        fragment_shader = shaders.fragment_shader_triangles, draw_type="triangles")
    xyz    = np.array([[0, 0.001, 0], [1, 1, 0], [-1, 1, 0], [-1, -1, 0], [1, -1, 0]], dtype=np.float32)*198
    colors = np.array([[1, 0, 0], [0, 1, 0], [ 0, 0, 1], [ 1,  1, 0], [1,  0, 1]], dtype=np.float32)
    radii = np.random.randint(20,50,5).astype(np.float32)
    data = GLContainer(xyz, colors, radii)
    test.load_data(data)
    wind = Gtk.Window()
    wind.add(test)

    wind.connect("delete-event", Gtk.main_quit)
    wind.connect("key-press-event", test.key_pressed)
    wind.connect("key-release-event", test.key_released)
    wind.show_all()
    Gtk.main()

