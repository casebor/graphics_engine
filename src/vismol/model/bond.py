#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np


class Bond:
    """ Class doc """
    
    def __init__(self, atom_i, atom_j, atom_index_i=None, atom_index_j=None, bond_order=1):
        """ Class initialiser """
        self.atom_i = atom_i
        self.atom_j = atom_j
        # - Remember that the "index" attribute refers to the numbering of atoms 
        # (it is not a zero base, it starts at 1 for the first atom)
        # these indices are zero base numbering (it starts at 0 for the first atom)
        if atom_index_i:
            self.atom_index_i = atom_index_i
        else:
            self.atom_index_i = atom_i.index-1
        
        if atom_index_j:
            self.atom_index_j = atom_index_j
        else:
            self.atom_index_j = atom_j.index-1
        
        self.bond_order = bond_order
        self.line_active = True
        self.stick_active = False
    
    def distance (self, frame=0):
        """ Function doc """
        vec = self.atom_i.coords(frame) - self.atom_j.coords(frame)
        return np.linalg.norm(vec)
