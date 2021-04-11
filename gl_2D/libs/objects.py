#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import numpy as np

class VMNode():
    """docstring for VMAtom"""
    
    def __init__ (self, nid, x, y, z, normal=None, color=None):
        self.id = nid
        self.x = x
        self.y = y
        self.z = z
        self.pos = np.array([x, y, z], dtype=float)
        if normal is not None:
            self.normal = np.array(normal, dtype=float)
        else:
            self.normal = self.pos / np.linalg.norm(self.pos)
        if color is None:
            self.color = np.array([0,1,0], dtype=float)
        else:
            self.color = np.array(color)
        # self.bonds = []
    
    def set_pos(self, new_pos):
        """ Function doc """
        assert len(new_pos) == 3
        self.x, self.y, self.z = new_pos
        self.pos = np.array(new_pos, dtype=float)
    
    def set_normal(self, new_normal):
        """ Function doc """
        assert len(new_normal) == 3
        # assert np.linalg.norm(new_normal) == 1.0
        self.normal = np.array(new_normal, dtype=float)
    
    # def add_bond(self, other_node):
    #     """ Function doc """
    #     if other_node not in self.bonds:
    #         self.bonds.append(other_node)
    
    # def del_bond(self, other_node):
    #     """ Function doc """
    #     if other_node in self.bonds:
    #         self.bonds.remove(other_node)

class VMAtom(VMNode):
    """docstring for VMAtom"""
    def __init__(self, gl_id, atom_id, x, y, z, name="UNK", radius=1.0, element="X",
                 parent=None, color=None):
        super().__init__(gl_id, x, y, z, color=color)
        self.atom_id = atom_id
        self.rad = radius
        self.name = name
        self.element = element
        self.parent = parent

class GLContainer():
    """docstring for GLContainer"""
    def __init__(self, obj_list: list, indexes=None):
        self.xyz = np.zeros((len(obj_list), 3), dtype=float)
        self.colors = np.zeros((len(obj_list), 3), dtype=float)
        self.radii = np.zeros((len(obj_list)), dtype=float)
        self.vel = 0.1
        for i, obj in enumerate(obj_list):
            self.xyz[i,:] = obj.pos
            self.colors[i,:] = obj.color
        if indexes is not None:
            self.indexes = np.array(indexes, dtype=np.uint32)
        else:
            self.indexes = None
    
    def build_dirs(self):
        self.dir = self.xyz / np.linalg.norm(self.xyz, axis=0)

    def update_pos(self, vel=0.1):
        to_add = self.dir * vel
        self.xyz += to_add
