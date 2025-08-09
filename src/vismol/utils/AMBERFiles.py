#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# TODO: Needs revision

import os
import numpy as np
from typing import List
from logging import getLogger
from vismol.core.vismol_object import VismolObject


logger = getLogger(__name__)


def load_netcdf4_file(vm_object: "VismolObject", infile: str) -> np.array:
    """
    Load frames from a NetCDF4 file.
    
    Args:
        vm_object (VismolObject): The VismolObject.
        infile (str): Path to the NetCDF4 file.
    
    Returns:
        np.ndarray: Array of frames.
    
    """
    from netCDF4 import Dataset
    f = Dataset(infile)
    size = len(vm_object.atoms) * 3
    frames = np.array(f.variables["coordinates"][:], dtype=np.float32)
    frames = frames.reshape((frames.shape[0], frames.shape[1] * frames.shape[2]))
    return frames


def load_amber_crd_file(vm_object: "VismolObject", infile: str) -> List[np.ndarray]:
    """
    Load frames from an AMBER CRD file.
    
    Args:
        vm_object (VismolObject): The VismolObject.
        infile (str): Path to the AMBER CRD file.
    
    Returns:
        List[np.ndarray]: List of frame arrays.
    
    """
    size = len(vm_object.atoms) * 3
    with open(infile, "r") as infile_handler:
        data = infile_handler.readlines()
        _line = data[1].split()
        if int(float(_line[0])) == size / 3:
            start = 2
        else:
            start = 1
        frames = []
        frame = []
        counter = 0
        for line in data[start:]:
            _line = line.split()
            for coord in _line:
                frame.append(float(coord))
                counter += 1
                if counter == size:
                    frame = np.array(frame, dtype=np.float32)
                    frames.append(frame)
                    frame = []
                    counter = 0
    return frames


def load_amber_topology_file(vm_session: "VismolSession", infile: str, gridsize: int = 3) -> "VismolObject":
    """
    Load AMBER topology information from a file.
    
    Args:
        vm_session (VismolSession): The VismolSession.
        infile (str): Path to the input file.
        gridsize (int, optional): Grid size (default is 3).
    
    Returns:
        None
    
    """
    with open(infile, "r") as infile_handler:
        text = infile_handler.read()
        text = text.split("%FLAG")
        atom_names = None
        atomic_numbers = None
        bonds_H = None
        bonds_noH = None
        residues_names = None
        pointers = None
        at_names = []
        at_numbers = []
        bonds = []
        bonds2 = []
        res_names = []
        res_indexes = []
        total_bonds = []
        
    for block1 in text:
        block2 = block1.split("\n")
        string = ""
        for block3  in block2[2:]:
            string += block3
        
        if "ATOM_NAME" in block1:
            atom_names = string
            for i in range(0, len(atom_names),4):
                at_name = atom_names[i:i+4]
                at_names.append(at_name.strip())
        
        elif "ATOMIC_NUMBER" in block1:
            atomic_numbers = string
            for i in range(0, len(atomic_numbers),8):
                atomic_number = atomic_numbers[i:i+8]
                at_numbers.append(int(atomic_number.strip()))
        
        elif "BONDS_INC_HYDROGEN" in block1:
            bonds_H = string
            for i in range(0, len(bonds_H),8):
                bond = bonds_H[i:i+8]
                bond = int(int(bond.strip())/3)
                bonds.append(bond)
            for i in range(0, len(bonds),3):
                total_bonds.append(bonds[i:i+2])
        
        elif "BONDS_WITHOUT_HYDROGEN" in block1:
            bonds_noH = string
            for i in range(0, len(bonds_noH),8):
                bond = bonds_noH[i:i+8]
                bond = int(int(bond.strip())/3)
                bonds2.append(bond)
            for i in range(0, len(bonds2),3):
                total_bonds.append(bonds2[i:i+2])
        
        elif "RESIDUE_LABEL" in block1:
            labels = string
            res_names_short = []
            for i in range(0, len(labels),4):
                label = labels[i:i+4]
                res_names_short.append(label.strip())
        
        elif "RESIDUE_POINTER" in block1:
            pointers = string.split()
            pointers2 = []
            for i in pointers:
                pointers2.append(int(i))
            pointers2.append(len(at_names)+1)
            idx = 0
            residue_counter = 1
            for atom in at_names:
                if idx < pointers2[residue_counter]-1:
                    res_indexes.append(residue_counter)
                    res_names.append(res_names_short[residue_counter-1])
                else:
                    residue_counter += 1
                    res_indexes.append(residue_counter)
                    res_names.append(res_names_short[residue_counter-1])
                idx += 1
    
    resi_counter  = 1
    atoms = []
    
    for index in range(len(at_names)):
        at_name  = at_names[index]
        at_pos = np.array([0.0,0.0,0.0])
        at_resn = res_names[index]
        at_resi = res_indexes[index]
        at_ch = "X"
        at_occup = 0.0
        at_bfactor = 0.0
        at_charge = 0.0
        at_symbol = None
        atoms.append({"index": index, "name": at_name, "resi": at_resi,
                      "resn": at_resn, "chain": at_ch, "symbol": at_symbol,
                      "occupancy": at_occup, "bfactor": at_bfactor,
                      "charge": at_charge})
    
    bonds_pair_of_indexes = []
    for bond in total_bonds:
        bonds_pair_of_indexes.append(bond[0])
        bonds_pair_of_indexes.append(bond[1])
    
    name = os.path.basename(filename)
    vm_object  = VismolObject(name=name, atoms=atoms, vm_session=vm_session,
                              bonds_pair_of_indexes=bonds_pair_of_indexes,
                              trajectory=[], auto_find_bonded_and_nonbonded=False)
    return vm_object

