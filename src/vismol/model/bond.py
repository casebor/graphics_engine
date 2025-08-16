#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np


class Bond:
    """ Class doc """
    
    def __init__(self, atom_i: "Atom", atom_j: "Atom", bond_order: np.int32=1):
        """ Class initialiser """
        self.atom_i = atom_i
        self.atom_j = atom_j
        # - Remember that the "index" attribute refers to the numbering of atoms 
        # (it is not a zero base, it starts at 1 for the first atom)
        # these indices are zero base numbering (it starts at 0 for the first atom)
        self.atom_index_i = atom_i.atom_id
        self.atom_index_j = atom_j.atom_id
        # if atom_index_i:
        # else:
        #     self.atom_index_i = atom_i.index-1
        
        # if atom_index_j:
        #     self.atom_index_j = atom_index_j
        # else:
        #     self.atom_index_j = atom_j.index-1
        
        self.bond_order = bond_order
        self.line_active = True
        self.stick_active = False
    
    def __hash__(self):
        return hash(frozenset({self.atom_i, self.atom_j}))
    
    def __eq__(self, other: "Bond") -> bool:
        return isinstance(other, Bond) and \
               frozenset({self.atom_i, self.atom_j}) == frozenset({other.atom_i, other.atom_j})
    
    def __repr__(self):
        return "Bond atoms: {}-{}\t IDs: {}-{}".format(self.atom_i, self.atom_j, self.atom_index_i, self.atom_index_j)
    
    def distance(self, frame: np.int32=0) -> np.float32:
        """ Function doc """
        vec = self.atom_i.get_coords_from_frame(frame) - self.atom_j.get_coords_from_frame(frame)
        return np.linalg.norm(vec)
