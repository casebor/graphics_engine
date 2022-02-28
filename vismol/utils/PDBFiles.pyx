#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
import os
import time
import multiprocessing
import numpy as np
cimport numpy as np
from logging import getLogger
from core.vismol_object import VismolObject
from model.molecular_properties import AtomTypes

logger = getLogger(__name__)


cpdef get_topology(infile):
    """ This function should return the list of tuples with the information of
        the atoms in the PDB file
    """
    initial = time.time()
    rawframe = ""
    with open(infile, "r") as pdb_file:
        for line in pdb_file:
            if line[:3] == "END":
                break
            elif (line[:4] == "ATOM") or (line[:6] == "HETATM"):
                rawframe += line
    atoms_info = get_list_of_atoms_from_rawframe(rawframe)
    logger.debug("Time used to parse the topology: {:>8.5f} secs".format(time.time() - initial))
    return atoms_info

cpdef get_coords(infile, atom_qtty, n_proc):
    """ Here we should focus on parsing the coordinates only
        This method of splitting the file will double the amount of RAM memory
        required. A for line in file could be used, however this approach will
        increase the time required :S
    """
    initial = time.time()
    raw_frames = []
    raw_frame = ""
    cdef int i = 0
    with open(infile, "r") as pdbf:
        for line in pdbf:
            if line[:3] == "END" and i != 0:
                raw_frames.append(raw_frame)
                raw_frame = ""
                i = 0
            elif (line[:4] == "ATOM") or (line[:6] == "HETATM"):
                raw_frame += line
                i += 1
    coords = get_coords_from_pdb_rawframes(raw_frames, len(raw_frames), atom_qtty, n_proc)
    logger.debug("Time used to parse the coordinates: {:>8.5f} secs".format(time.time() - initial))
    return coords

cpdef _get_ranges(traj_size, batch_size):
    cdef int i = 0
    ranges = []
    while (i + batch_size) < traj_size:
        ranges.append((i, i + batch_size))
        i += batch_size
    ranges.append((i, traj_size))
    return ranges

cpdef get_coords_from_pdb_rawframes(rawframes, traj_size, atom_qtty, n_proc):
    """ Function doc """
    with multiprocessing.Pool(processes=n_proc) as pool:
        processing = []
        batch_size = np.int32(np.ceil(traj_size / n_proc))
        for r in _get_ranges(traj_size, batch_size):
            processing.append(pool.apply_async(get_pdb_frame_coords,
                              args=(rawframes[r[0]:r[1]], atom_qtty)))
        frames = [p.get() for p in processing]
    coords = np.empty([0, atom_qtty, 3], dtype=np.float32)
    for f in frames:
        coords = np.vstack((coords, f))
    return coords

cpdef get_pdb_frame_coords(frames, atom_qtty):
    """ Function doc """
    cdef int i = 0
    cdef int j = 0
    coords = np.empty([len(frames), atom_qtty, 3], dtype=np.float32)
    for frame in frames:
        j = 0
        lines = frame.strip().split("\n")
        for line in lines:
            x = np.float32(line[30:38])
            y = np.float32(line[38:46])
            z = np.float32(line[46:54])
            coords[i,j,:] = x, y, z
            j += 1
        i += 1
    return coords


cpdef load_pdb_file(infile, vismol_session, gridsize=3, frames_only=False):
    """ The longest covalent bond I can find is the bismuth-iodine single bond.
        The order of bond lengths is single > double > triple.
        The largest atoms should form the longest covalent bonds. 
        So we look at atoms in the lower right corner of the Periodic Table.
        The most likely candidates are Pb, Bi, and I.
        The experimental bond lengths are:
        Bi-I = 281 pm; Pb-I = 279 pm; I-I = 266.5 pm.
        So the polar covalent Bi-I bond is the longest covalent measured so far.
    """
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
            atoms_info = get_list_of_atoms_from_rawframe(rawframes[i])
            i += 1
    
    name = os.path.basename(infile)
    vismol_object = VismolObject(vismol_session, atoms_info, name=name, trajectory=frames)
    return vismol_object

cpdef get_list_of_frames_from_pdb_rawframes(rawframes=None):
    """ Function doc """
    n_processor = multiprocessing.cpu_count()
    pool        = multiprocessing.Pool(n_processor)
    frames      = pool.map(get_pdb_frame_coordinates, rawframes)
    framesout   = []
    
    for frame in frames:
        if frame:
            frame = np.array(frame, dtype=np.float32)
            framesout.append(frame)
    
    return framesout

cpdef get_list_of_atoms_from_rawframe(rawframe):
    """ Function doc """
    pdb_file_lines = rawframe.split("\n")
    atoms = []
    cdef int index = 1
    for line in pdb_file_lines:
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            at_index  = np.int32(line[7:11])
            at_name   = line[12:16].strip()
            at_resn   = line[17:20].strip()
            at_ch     = line[20:22].strip()
            at_resi   = np.int32(line[22:27])
            at_occup   = np.float32(line[54:60])
            at_bfactor = np.float32(line[60:66])
            at_charge  = 0.0
            if at_index != index:
                logger.debug("Atom index {} is different from number of atoms in rawframe".format(at_index))
            atoms.append({"name"     : at_name,
                          "resn"     : at_resn,
                          "chain"    : at_ch,
                          "resi"     : at_resi,
                          "occupancy": at_occup,
                          "bfactor"  : at_bfactor,
                          "charge"   : at_charge,
                          "index"    : index})
            index += 1
    return atoms

cpdef get_pdb_frame_coordinates(str frame):
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
