#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
import os
import time
import multiprocessing
import numpy as np
cimport numpy as np
from core.vismol_object import VismolObject
from model.molecular_properties import AtomTypes


cdef load_pdb_file(infile, vismol_session, gridsize=3, frames_only=False):
    """ The longest covalent bond I can find is the bismuth-iodine single bond.
        The order of bond lengths is single > double > triple.
        The largest atoms should form the longest covalent bonds. 
        So we look at atoms in the lower right corner of the Periodic Table.
        The most likely candidates are Pb, Bi, and I.
        The experimental bond lengths are:
        Bi-I = 281 pm; Pb-I = 279 pm; I-I = 266.5 pm.
        So the polar covalent Bi-I bond is the longest covalent measured so far.
    """
    #-------------------------------------------------------------------------------------------
    #                                P D B     P A R S E R 
    #-------------------------------------------------------------------------------------------
    initial = time.time()
    with open(infile, "r") as pdb_file:
        pdbtext = pdb_file.read()
        if "ENDMDL" in pdbtext:
            rawframes = pdbtext.split("ENDMDL")
        else:
            rawframes = pdbtext.split("END")
        frames = get_list_of_frames_from_pdb_rawframes(rawframes)
        
        if frames_only:
            return frames
        i = 0
        atoms_info = []
        # TODO: NOT SURE IF THIS IS A CORRECT USE OF WHILE
        while atoms_info == []:
            atoms_info = get_list_of_atoms_from_rawframe(rawframes[i], gridsize)
            i += 1
    
    name = os.path.basename(infile)
    final = time.time()
    # print ("P D B     P A R S E R end -  total time: ", final - initial, "\n")
    
    vismol_object = VismolObject(vismol_session, atoms_info, name=name, trajectory=frames)
    return vismol_object

cdef get_list_of_atoms_from_rawframe(rawframe, gridsize=3):
    """ Function doc 
    """
    at = AtomTypes()
    pdb_file_lines = rawframe.split("\n")
    atoms = []
    
    cdef int index = 0
    for line in pdb_file_lines:
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            at_name   = line[12:16].strip()
            at_pos    = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            at_resi   = int(line[22:27])
            at_resn   = line[17:20].strip()
            at_ch     = line[20:22].strip()
            at_occup   = float(line[54:60])
            at_bfactor = float(line[60:66])
            at_charge  = 0.0
            at_symbol = line[70:].strip()
           
            if at_symbol not in at.ATOM_TYPES.keys():
                at_symbol = None #at.get_symbol(at_name)
            
            atoms.append({"index"    : index,
                          "name"     : at_name,
                          "resi"     : at_resi,
                          "resn"     : at_resn,
                          "chain"    : at_ch,
                          "symbol"   : at_symbol,
                          "occupancy": at_occup,
                          "bfactor"  : at_bfactor,
                          "charge"   : at_charge})
            index += 1
    return atoms

cdef get_list_of_frames_from_pdb_rawframes(rawframes=None):
    """ Function doc
    """
    n_processor = multiprocessing.cpu_count()
    pool        = multiprocessing.Pool(n_processor)
    frames      = pool.map(get_pdb_frame_coordinates, rawframes)
    framesout   = []
    
    for frame in frames:
        if frame:
            frame = np.array(frame, dtype=np.float32)
            framesout.append(frame)
    
    return framesout

cdef get_pdb_frame_coordinates(str frame):
    """ Function doc """
    pdb_file_lines = frame.split("\n")
    frame_coordinates = []
    
    for line in pdb_file_lines:
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            frame_coordinates.append(float(line[30:38]))
            frame_coordinates.append(float(line[38:46]))
            frame_coordinates.append(float(line[46:54]))
    
    if len(frame_coordinates) == 0:
        return False
    else:
        return frame_coordinates
