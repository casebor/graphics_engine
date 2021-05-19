#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import numpy as np

class GLContainer():
    """docstring for GLContainer"""
    def __init__(self, coords:np.ndarray, colors:np.ndarray, radii:np.ndarray,
                 direcs=None, velocity=1.0):
        self.xyz = coords
        self.colors = colors
        self.radii = radii
        self.velocity = velocity
        if direcs is not None:
            self.direc = np.array(direcs, dtype=float)
        else:
            self.direc = (self.xyz.T / np.linalg.norm(self.xyz, axis=1)).T
        self.borders = np.array([400.,400.])
    
    def update_pos(self, vel=0.1):
        to_add = self.direc[:,:2] * vel
        # flags = self.xyz[:,:2] + to_add > self.borders
        # flags = flags.any(axis=1)
        # for i,flag in enumerate(flags):
        #     if flag:
        #         if self.xyz[i,0] + to_add[i,0] > self.borders[0]:
        #             self.direc[i,0] *= -1
        #         if self.xyz[i,1] + to_add[i,1] > self.borders[1]:
        #             self.direc[i,1] *= -1
        # self.xyz[:,:2] += to_add
        
        flags = self.xyz[:,:2] + to_add > self.borders
        flags = flags.astype(int)
        self.xyz[:,:2] += to_add*flags
        flags = self.xyz[:,:2] + to_add < np.zeros(2,dtype=float)
        flags = flags.astype(int)*-1
        self.xyz[:,:2] += to_add*flags
