#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import numpy as np

class Node2D():
    """docstring for Node2D"""
    def __init__(self, x=0.0, y=0.0, radius=1.0, color=[0.0, 0.0, 0.0]):
        self.x = x
        self.y = y
        self.rad = radius
        self.color = np.array(color, dtype=np.float64)
        self.position = np.array([x, y], dtype=np.float64)
        self.velocity = None
        self.direction = None
    
    def set_velocity(self, velocity):
        """ Docs for set_velocity """
        assert(isinstance(velocity, float))
        self.velocity = velocity
    
    def set_direction(self, direction):
        """ Docs for set_direction """
        direction = np.array(direction, dtype=np.float64)
        assert direction.shape[0] == 2
        self.direction = direction / np.linalg.norm(direction)
    
    def set_position(self, position):
        """ Docs for set_position """
        position = np.array(position, dtype=np.float64)
        assert position.shape[0] == 2
        self.position = position
        self.x = position[0]
        self.y = position[1]
    
    def update_position(self):
        """ Docs for update_position """
        assert self.velocity is not None
        assert self.direction is not None
        self.position += self.velocity * self.direction
        self.x = self.position[0]
        self.y = self.position[1]

class Box2D():
    """docstring for Box2D"""
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.nodes = []
        self.dt = 0.1

    def add_node(self, node):
        """ Docs for add_node """
        assert isinstance(node, Node2D)
        if node not in self.nodes:
            self.nodes.append(node)

    def set_width(self, width):
        """ Docs for set_width """
        self.width = width

    def set_height(self, height):
        """ Docs for set_height """
        self.height = height
