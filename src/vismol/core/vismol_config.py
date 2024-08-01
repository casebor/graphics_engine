#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import json


class VismolConfig:
    """
    Class containing the configuration options for Vismol.
    TODO: Some configuration options are not implemented yet, but they are used
    in the code, need to update all the parameters :(
    
    Attributes:
        vm_session (VismolSession): The Vismol session for which this
                configuration is applied.
        gl_parameters (Dict): Dictionary of parameters for the OpenGL engine.
        n_procs (int): Number of cores to use.
        represetations_available (List): List of represetations available for
                the session.
    """
    
    def __init__ (self, vm_session: "VismolSession"):
        """
        Args:
            vm_session (VismolSession): The VismolSession object.
        """
        self.vm_session = vm_session
        self.gl_parameters = {"antialias": True,
                              "background_color": [0.0, 0.0, 0.0, 1.0],
                              "bond_tolerance": 1.4,
                              "center_on_coord_sleep_time": 0.01,
                              "color_type": 0,
                              "dashed_dist_lines_color": [0.4, 0.4, 0.4, 1.0],
                              "dot_sel_size": 1.5,
                              "dot_size": 2,
                              "dot_type": 1,
                              "dots_size": 2,
                              "impostor_type": 0,
                              "field_of_view": 10,
                              "gridsize": 0.8,
                              "light_ambient_coef": 0.4,
                              "light_color": [ 1.0, 1.0, 1.0, 1.0],
                              "light_intensity": [0.6, 0.6, 0.6],
                              "light_position": [0, 0, 10.0],
                              "light_shininess": 5.5,
                              "light_specular_color": [1.0, 1.0, 1.0],
                              "line_color": 0,
                              "line_width": 3,
                              "line_width_selection": 80,
                              "line_type": 0,
                              "maxbond": 2.4,
                              "picking_dots_color": [0.0, 1.0, 1.0],
                              "picking_dots_safe": True,
                              "pk_label_color": [1.0, 1.0, 1.0, 1.0],
                              "pk_dist_label_color": [1.0, 1.0, 0.0, 1.0],
                              "ribbon_color": 0,
                              "ribbon_width": 1000,
                              "ribbon_width_selection": 100,
                              "ribbon_type": 2,
                              "scroll_step": 0.9,
                              "sphere_quality": 2,
                              "sphere_scale": 0.20,
                              "sphere_type": 0,
                              "sticks_color": 0,
                              "sticks_radius": 0.10,
                              "sticks_type": 0}
        self.n_proc = 2
        self.representations_available = {"dots", "lines", "nonbonded",
            "impostor","dash", "sticks", "spheres", "ribbons", "dynamic",
            "vdw_spheres", "picking_spheres", "static_freetype"}
    
    
    def save_easyhybrid_config(self) -> None:
        """
        Save the configuration to a JSON file.
        
        """
        os.makedirs(os.path.join(os.environ["HOME"], ".VisMol"), exist_ok=True)
        config_path = os.path.join(os.environ["HOME"], ".VisMol", "VismolConfig.json")
        with open(config_path, "w") as config_file:
            json.dump(self.gl_parameters, config_file, indent=2)
    
    
    def load_easyhybrid_config(self, config_path: str) -> None:
        """
        Load the configuration from a JSON file.
        
        Args:
            config_path (str): Path to the configuration file.
        """
        if not os.path.isfile(config_path):
          config_path = os.path.join(os.environ["HOME"], ".VisMol", "VismolConfig.json")
        with open(config_path, "r") as config_file:
            self.gl_parameters = json.load(config_file)
