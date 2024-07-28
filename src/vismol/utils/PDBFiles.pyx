#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import time
import multiprocessing
import numpy as np
cimport numpy as np
from logging import getLogger
from vismol.model.atom import Atom
from vismol.model.residue import Residue
from vismol.model.chain import Chain
from vismol.model.molecule import Molecule
from vismol.core.vismol_object import VismolObject


logger = getLogger(__name__)


cpdef get_list_of_atoms_from_rawframe(raw_frames):
    """ Function doc """
    pdb_file_lines = raw_frames.split("\n")
    atoms = []
    cdef int index = 1
    for line in pdb_file_lines:
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            at_index = index
            at_name = line[12:16].strip()
            at_resn = line[17:20].strip()
            at_ch = line[20:22].strip()
            at_resi = np.int32(line[22:27])
            at_occup = np.float32(line[54:60])
            at_bfactor = np.float32(line[60:66])
            at_charge = 0.0
            atoms.append({"name": at_name, "resn": at_resn, "chain": at_ch,
                          "resi": at_resi, "occupancy": at_occup,
                          "bfactor": at_bfactor, "charge": at_charge,
                          "index": index})
            index += 1
    return atoms


cpdef get_topology(infile):
    """ This function should return the list of tuples with the information of
        the atoms in the PDB file
    """
    initial = time.time()
    with open(infile, "r") as pdb_file:
        pdbtext = pdb_file.read()
        if "ENDMDL" in pdbtext:
            raw_frames = pdbtext.split("ENDMDL")
        else:
            raw_frames = pdbtext.split("END")
    atoms_info = get_list_of_atoms_from_rawframe(raw_frames[0])
    logger.debug("Time used to parse the topology: {:>8.5f} secs".format(time.time() - initial))
    return atoms_info, raw_frames


cpdef _get_frame_coords(frames, atom_qtty):
    """ Function doc """
    cdef int i = 0
    cdef int j = 0
    coords = np.empty([len(frames), atom_qtty, 3], dtype=np.float32)
    for frame in frames:
        j = 0
        lines = frame.strip().split("\n")
        for line in lines:
            if line[:4] == "ATOM" or line[:6] == "HETATM":
                x = np.float32(line[30:38])
                y = np.float32(line[38:46])
                z = np.float32(line[46:54])
                coords[i,j,:] = x, y, z
                j += 1
        i += 1
    return coords


cpdef _get_ranges(traj_size, batch_size):
    cdef int i = 0
    ranges = []
    while (i + batch_size) < traj_size:
        ranges.append((i, i + batch_size))
        i += batch_size
    ranges.append((i, traj_size))
    return ranges


cpdef _get_coords_parallel(raw_frames, traj_size, atom_qtty, n_proc):
    """ Function doc """
    with multiprocessing.Pool(processes=n_proc) as pool:
        processing = []
        batch_size = np.int32(np.ceil(traj_size / n_proc))
        for r in _get_ranges(traj_size, batch_size):
            processing.append(pool.apply_async(_get_frame_coords,
                              args=(raw_frames[r[0]:r[1]], atom_qtty)))
        frames = [p.get() for p in processing]
    coords = np.empty([0, atom_qtty, 3], dtype=np.float32)
    for f in frames:
        coords = np.vstack((coords, f))
    return coords


cpdef get_coords_from_raw_frames(raw_frames, atom_qtty, n_proc):
    """ Here we should focus on parsing the coordinates only
        This method of splitting the file will double the amount of RAM memory
        required. A for line in file could be used, however this approach will
        increase the time required :S
    """
    initial = time.time()
    if not "ATOM" in raw_frames[0] and not "HETATM" in raw_frames[0]:
        raw_frames.pop(0)
    if not "ATOM" in raw_frames[-1] and not "HETATM" in raw_frames[-1]:
        raw_frames.pop(-1)
    coords = _get_coords_parallel(raw_frames, len(raw_frames), atom_qtty, n_proc)
    logger.debug("Time used to parse the coordinates: {:>8.5f} secs".format(time.time() - initial))
    return coords


cpdef get_coords_from_file(infile, atom_qtty, n_proc):
    """ Here we should focus on parsing the coordinates only
        This method of splitting the file will double the amount of RAM memory
        required. A for line in file could be used, however this approach will
        increase the time required :S
    """
    initial = time.time()
    cdef int i = 0
    with open(infile, "r") as pdbf:
        pdbtext = pdbf.read()
        if "ENDMDL" in pdbtext:
            raw_frames = pdbtext.split("ENDMDL")
        else:
            raw_frames = pdbtext.split("END")
    if not "ATOM" in raw_frames[0] and not "HETATM" in raw_frames[0]:
        raw_frames.pop(0)
    if not "ATOM" in raw_frames[-1] and not "HETATM" in raw_frames[-1]:
        raw_frames.pop(-1)
    coords = _get_coords_parallel(raw_frames, len(raw_frames), atom_qtty, n_proc)
    logger.debug("Time used to parse the coordinates: {:>8.5f} secs".format(time.time() - initial))
    return coords


cpdef load_pdb_file(vm_session, infile):
    vm_object_name = os.path.basename(infile)
    topology, rawframes = get_topology(infile)
    vm_object = VismolObject(vm_session, len(vm_session.vm_objects_dic),
                             name=vm_object_name)
    unique_id = len(vm_session.atom_dic_id)
    initial = time.time()
    atom_id = 0
    
    for _atom in topology:
        if vm_object.molecule is None:
            vm_object.molecule = Molecule(vm_object)
        _mol = vm_object.molecule
        
        if _atom["chain"] not in _mol.chains:
            _mol.chains[_atom["chain"]] = Chain(vm_object, name=_atom["chain"],
                                                molecule=_mol)
        _chain = _mol.chains[_atom["chain"]]
        
        if _atom["resi"] not in _chain.residues:
            _res = Residue(vm_object, name=_atom["resn"], index=_atom["resi"],
                           chain=_chain)
            _chain.residues[_atom["resi"]] = _res
        _residue = _chain.residues[_atom["resi"]]
        
        atom = Atom(vm_object, name=_atom["name"], index=_atom["index"],
                    residue=_residue, chain=_chain, atom_id=atom_id,
                    occupancy=_atom["occupancy"], bfactor=_atom["bfactor"],
                    charge=_atom["charge"], molecule=_mol)
        atom.unique_id = unique_id
        atom._generate_atom_unique_color_id()
        _residue.atoms[atom_id] = atom
        _mol.atoms[atom_id] = atom
        atom_id += 1
        unique_id += 1
    
    logger.debug("Time used to build the tree: {:>8.5f} secs".format(time.time() - initial))
    _mol.frames = get_coords_from_raw_frames(rawframes, atom_id,
                                             vm_session.vm_config.n_proc)
    _mol.geom_center = np.mean(_mol.frames[0], axis=0)
    _mol.find_bonded_and_nonbonded_atoms()
    return vm_object


