#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_config.py
#  
#  Copyright 2022 Fernando <fernando@winter>
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
import json

class VismolConfig:
    """ Class doc """
    
    def __init__ (self, vismol_session):
        """ Class initialiser """
        self.vismol_session = vismol_session
        self.gl_parameters = {"background_color": [0.0, 0.0, 0.0, 1.0],
                              "color_type": 0,
                              "dot_size": 20,
                              "dots_size": 2,
                              "dot_type": 1,
                              "dot_sel_size": 1.5,
                              "line_width": 2,
                              "line_width_selection": 80,
                              "line_type": 0,
                              "line_color": 0,
                              "ribbon_width": 1000,
                              "ribbon_width_selection": 100,
                              "ribbon_type": 1,
                              "ribbon_color": 0,
                              "sphere_type": 0,
                              "sphere_scale": 1.0,
                              "sphere_quality": 1,
                              "impostor_type": 1,
                              "stick_radius": 3.5,
                              "stick_color": 0,
                              "stick_type": 0,
                              "antialias": True,
                              "scroll_step": 0.9,
                              "field_of_view": 10,
                              "light_position": [-2.5, -2.5, 3.0],
                              "light_color": [ 1.0, 1.0, 1.0,1.0],
                              "light_ambient_coef": 0.4,
                              "light_shininess": 5.5,
                              "light_intensity": [0.6, 0.6, 0.6],
                              "light_specular_color": [1.0, 1.0, 1.0],
                              "center_on_coord_sleep_time": 0.001,
                              "gridsize": 0.8,
                              "maxbond": 2.4,
                              "bond_tolerance": 1.4,
                              "picking_dots_color": [0.0, 1.0, 1.0]}
        self.n_proc = 2
        # self.representations_available = {"dots", "lines", "nonbonded", "dotted_lines",
        #                                   "ribbon", "sticks", "spheres", "impostor",
        #                                   "surface", "cartoon", "freetype",
        #                                   "picking_dots"}
        self.representations_available = {"dots", "lines", "nonbonded", "picking_dots"}
    
    
    def save_easyhybrid_config(self):
        """ Function doc """
        config = os.path.join(os.environ["HOME"], ".VisMol", "VismolConfig.json")
        with open(config, "w") as config_file:
            json.dump(self.gl_parameters, config_file, indent=2)
    
    def load_easyhybrid_config(self):
        """ Function doc """
        config = os.path.join(os.environ["HOME"], ".VisMol", "VismolConfig.json")
        with open(config, "r") as config_file:
            self.gl_parameters = json.load(config_file)
    
