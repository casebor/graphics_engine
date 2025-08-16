#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class Chain:
    """ Class doc """
    
    def __init__ (self, vismol_object, name=None, molecule=None, residues=None, label=None):
        """ Class initialiser """
        self.vm_object = vismol_object
        self.molecule = molecule
        if residues is None:
            self.residues = {}
        else:
            self.residues = residues
        self.residues_by_index = {}
        self.backbone = []
        self.name = name
        self.backbone_pair_indexes_full = []
        self.backbone_pair_indexes_show = []
        self.has_altloc = False
    
    def get_CA_list(self):
        """ Function doc """
        self.backbone = []
        
        for residue in self.residues:
            if residue.isProtein:
                for atom in residue.atoms:
                    if atom.name == "CA":
                        self.backbone.append(atom)
                    else:
                        pass
        
        return self.backbone
    
    def get_secundary_structure(_):
        """ Function doc """
        return None
    
    def return_name (self):
        """ Function doc """
        return self.name
