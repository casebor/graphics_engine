#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from vismol.model.molecular_properties import solvent_dictionary
from vismol.model.molecular_properties import residues_dictionary


class Residue:
    """ Class doc """
    
    def __init__(self, vismol_object, name="UNK", index=None, atoms=None, chain=None):
        """ Class initialiser """
        self.vm_object = vismol_object
        self.name = name
        self.index = index
        if atoms is None:
            self.atoms = {}
        else:
            self.atoms = atoms
        self.chain = chain
        self.is_protein = False
        self.is_solvent = False
        self._is_protein()
        self.topology = {}
        self.has_altloc = False
    
    def _is_protein(self):
        """ Function doc """
        # is it a salvent molecule?
        if self.name in solvent_dictionary.keys():
            self.is_solvent = True
        
        
        else: # is it a protein residue?
            if self.name in residues_dictionary.keys():
                self.is_protein = True
            else:
                # . If the residue code is not in the dictionary, 
                #   then check whether the N, CA and C atoms are present
                N  = False
                CA = False
                C  = False
                
                for index, atom in self.atoms.items():
                    #print (atom.name)
                    if atom.name == "N":
                        N  = True
                    if atom.name == "C":
                        CA = True
                    if atom.name == "C":
                        C  = True
                
                #print('N,CA,C',N,CA, C , self.name)
                if N == True and CA == True and C == True:
                    self.is_protein = True
                    #print('True - N,CA,C',N,CA,C, self.name )
            

    
    def geometry_center(self, frame=0):
        """ Function doc """
        if frame > len(self.vm_object.frames)-1:
            frame = len(self.vm_object.frames)-1
        gc = np.zeros(3, dtype=np.float32)
        for atom in self.atoms.values():
            gc += atom.get_coords_from_frame(frame)
        gc /= len(self.atoms.values())
        return gc
    
    def get_center_of_mass(self, mass=False, frame=0):
        """ Function doc """
        frame_size = len(self.vm_object.frames)-1
        
        if frame <= frame_size:
            pass
        else:
            frame = frame_size
        
        total = len(self.atoms)
        sum_x = 0.0
        sum_y = 0.0
        sum_z = 0.0
        
        for atom in self.atoms.values():
            coord = atom.get_coords_from_frame(frame)
            sum_x += coord[0]
            sum_y += coord[1]
            sum_z += coord[2]
        
        self.geom_center = np.array([sum_x / total,
                                     sum_y / total, 
                                     sum_z / total])
    
    def get_phi_and_psi(self):
        """ Function doc """
        if self.is_protein:
            dihedral_atoms = { 
                             } 
            print(self.name,self.index) 
            for atom in self.atoms:
                print(self.name, atom.name, atom.symbol, atom.get_coords_from_frame(), atom.bonds, atom.connected2 )
        else:
            pass
