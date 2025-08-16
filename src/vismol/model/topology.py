#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np


class Topology():
    """ Class doc """
    
    def __init__(self, molecule: "Molecule"=None, single_bonds: list=None,
                 double_bonds: list=None, triple_bonds: list=None,
                 arom_rings: list=None, arom_bonds: list=None):
        self.molecule = molecule
        self.single_bonds = single_bonds
        self.double_bonds = double_bonds
        self.triple_bonds = triple_bonds
        self.arom_rings = arom_rings
        self.arom_bonds = arom_bonds
        self.c_alpha_bonds = None
        # TODO: Temporary attribute
        self.bonds_pair_of_indexes = None
        self.non_bonded_atoms = None

