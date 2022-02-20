#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import os
import time
import multiprocessing
import numpy as np
from core.vismol_object import VismolObject


def load_xyz_file(infile=None, vismol_session=None, gridsize=3):
    """ Function doc """
    print ("\nstarting: parse_mol2")
    #at  =  MolecularProperties.AtomTypes()

    #initial = time.time()

    with open(infile, "r") as xyz_file:
        
        #pdbtext          = xyz_file.readlines()
        #totalsize        = len(pdbtext)
        #framesize        = int(pdbtext[0])
        #number_of_frames = totalsize/(framesize+2)
        
        
        xyz_lines        = xyz_file.readlines()
        total_size       = len(xyz_lines)    
        xyz_model_size   = int(xyz_lines[0])
        model_size       = xyz_model_size+2
        number_of_frames = total_size/(xyz_model_size+2)

        
        frame0           = xyz_lines[2: model_size]
        
        
        atoms, frames    = get_atom_list_from_xyz_frame (raw_atoms= frame0, frame = True, gridsize = 3)#, at = at)
        
        models     = []
        
        for i in range(0, int(number_of_frames)): 
            #print (i*model_size) , (i+1)*model_size
            model = xyz_lines[(i*model_size) : (i+1)*model_size]
            models.append(model)
        
        n = 1
        for model_i  in models[1:]:
            frame = []
            for line in model_i[2:]:
                line2 = line.split()
                #print line2, [float(line[1]), float(line[2]), float(line[3])]
                #at_name  = line2[0].strip()
                #at_pos   = np.array([float(line2[1]), float(line2[2]), float(line2[3])])
                frame.append(float(line2[1]))
                frame.append(float(line2[2]))
                frame.append(float(line2[3]))
                
            frame = np.array(frame, dtype=np.float32)
            frames.append(frame)
            n += 1

    
    #-------------------------------------------------------------------------------------------
    #                                Bonded and NB lists 
    #-------------------------------------------------------------------------------------------
    #atoms, bonds_full_indexes, bonds_pair_of_indexes, NB_indexes_list = cdist.generete_full_NB_and_Bonded_lists(atoms)
    #-------------------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------------------
    #                         Building   V I S M O L    O B J
    #-------------------------------------------------------------------------------------------
    name = os.path.basename(infile)
    #print (atoms)
    #print (frames)
    vismol_object  = VismolObject.VismolObject(name        = name, 
                                               atoms       = atoms, 
                                               vismol_session   = vismol_session, 
                                               trajectory  = frames)
    
    
    #vismol_object._generate_atomtree_structure()
    #vismol_object._generate_atom_unique_color_id()
    #vismol_object.index_bonds       = bonds_full_indexes
    #vismol_object.index_bonds_pairs = bonds_pair_of_indexes
    #vismol_object.non_bonded_atoms  = NB_indexes_list
    #-------------------------------------------------------------------------------------------
    return vismol_object

def get_atom_list_from_xyz_frame(raw_atoms, frame = True, gridsize = 3):#, at = None):
    """ Function doc """
    #nCPUs =  multiprocessing.cpu_count()
    #pool  = multiprocessing.Pool(nCPUs)
    #pdb_file_lines  = frame.split("\n")   
    #atoms = (pool.map(parse_pdb_line, pdb_file_lines))
    
    atoms  = []
    frames = []
    frame_coordinates = []
    #print (raw_atoms)
    index  = 0
    for line in raw_atoms:
        line = line.split()
        if len(line) > 1:
            #print (line) 
            #index    = int(line[0])-1
            
            at_name = line[0]
            
            at_pos  = np.array([float(line[1]), float(line[2]), float(line[3])])
            
            at_resi = 1
            
            at_resn = "UNK"
            
            at_ch   = "X"          
            
            at_occup   = 0.0     #occupancy
            at_bfactor = 0.0
            at_charge  = 0.0
            

            #at_symbol = line[5].split(".")
            at_symbol = None
            #cov_rad   = at.get_cov_rad (at_name)
            #gridpos  = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
            
            #atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])

            #atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
            
            #print (index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge )
            atoms.append({
                          "index"      : index      , 
                          "name"       : at_name    , 
                          "resi"       : at_resi    , 
                          "resn"       : at_resn    , 
                          "chain"      : at_ch      , 
                          "symbol"     : at_symbol  , 
                          "occupancy"  : at_occup   , 
                          "bfactor"    : at_bfactor , 
                          "charge"     : at_charge   
                          })
            
            
            index += 1

            frame_coordinates.append(float(line[1]))
            frame_coordinates.append(float(line[2]))
            frame_coordinates.append(float(line[3]))
    
    frame_coordinates = np.array(frame_coordinates, dtype=np.float32)
    frames.append(frame_coordinates)
    
    #print (frames)
    #print (atoms)
    return atoms, frames#, coords

def get_bonds(raw_bonds):
    """ Function doc """
    index_bonds              = []
    index_bonds_pairs        = []
    index_bonds_pairs_orders = []
    
    #print (raw_bonds)
    #print ("Obtain bonds from original MOL2 file")
    for line in raw_atoms:
        line = line.split()
        if len(line) == 4:
            index    = int(line[0])            
            atom1    = int(line[1]-1)
            atom2    = int(line[2]-1)
            order    = line[3]

            index_bonds      .append(atom1)
            index_bonds      .append(atom2)
            index_bonds_pairs.append([atom1,atom2])
            
            index_bonds_pairs_orders.append(order)

    return [index_bonds, index_bonds_pairs]

