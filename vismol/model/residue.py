#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import numpy as np
import model.Vectors as vectors
from model.molecular_properties import solvent_dictionary
from model.molecular_properties import residues_dictionary


class Residue:
    """ Class doc """
    
    def __init__ (self, atoms   = None, 
                        name    = 'UNK', 
                        index   = None,
                        chain   = None,
                        vismol_object = None):
        """ Class initialiser """
        self.atoms     = []
        self.resi      = index
        self.resn      = name
        self.chain     = chain
        self.vismol_object   = vismol_object
        self.isProtein = False
        self.isSolvent = False
        self.mass_center = None
        self.is_protein ()


        self.topology = {
                         
                        }
        
        
    def get_center_of_mass (self, mass = False, frame = 0):
        """ Function doc """
        
        frame_size = len(self.vismol_object.frames)-1
        
        if frame <= frame_size:
            pass
        else:
            frame = frame_size
        
        total = len(self.atoms)
        
        #coord = [0,0,0]
        sum_x = 0.0
        sum_y = 0.0
        sum_z = 0.0
        
        for atom in self.atoms:
            coord = atom.coords (frame)
            sum_x += coord[0]
            sum_y += coord[1]
            sum_z += coord[2]
        
        self.mass_center = np.array([sum_x / total,
                                     sum_y / total, 
                                     sum_z / total])


    def is_protein (self):
        """ Function doc """
        #residues_dictionary = MolecularProperties.#self.vismol_object.vismolSession.vConfig.residues_dictionary
        #solvent_dictionary  = MolecularProperties.#self.vismol_object.vismolSession.vConfig.solvent_dictionary
        # is it a protein residue?
        if self.resn in residues_dictionary.keys():
            self.isProtein = True
        else:
            self.isProtein = False
        
        # is it a salvent molecule?
        if self.resn in solvent_dictionary.keys():
            self.isSolvent = True
        else:
            self.isSolvent = False
        
        #return self.isProtein
    
    def get_phi_and_psi (self):
        """ Function doc """
        if self.isProtein:
            
            dihedral_atoms = { 
                              
            
                             } 
            print(self.resn,self.resi) 
            for atom in self.atoms:
                print(self.resn, atom.name, atom.symbol, atom.coords(), atom.bonds, atom.connected2 )
        
        else:
            pass
            