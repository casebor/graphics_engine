#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  parser.py
#  
#  Copyright 2022 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import os
import time
import numpy as np
from logging import getLogger
from utils import AUXFiles
from utils import GROFiles
from utils import MOL2Files
from utils import PDBFiles
from utils import PSFFiles
from utils import AMBERFiles
from utils import XYZFiles
from model.atom import Atom
from model.residue import Residue
from model.chain import Chain
from core.vismol_object import VismolObject

logger = getLogger(__name__)

def _load_aux_file(vismol_session, infile):
    """ Function doc """
    vismol_object = AUXFiles.load_aux_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_gro_file(vismol_session, infile):
    """ Function doc
    """
    vismol_object = GROFiles.load_gro_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_mol2_file(vismol_session, infile):
    """ Function doc """
    vismol_object = MOL2Files.load_mol2_files(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_pdb_file(vismol_session, infile):
    """ Function doc
    """
    vm_object_name = os.path.basename(infile)
    topo, rawframes = PDBFiles.get_topology(infile)
    vm_object = VismolObject(vismol_session, len(vismol_session.vm_objects_dic), name=vm_object_name)
    vm_object.set_model_matrix(vismol_session.vm_glcore.model_mat)
    unique_id = len(vismol_session.atom_dic_id)
    initial = time.time()
    atom_id = 0
    for _atom in topo:
        if _atom["chain"] not in vm_object.chains.keys():
            vm_object.chains[_atom["chain"]] = Chain(vm_object, name=_atom["chain"])
        _chain = vm_object.chains[_atom["chain"]]
        
        if _atom["resi"] not in _chain.residues.keys():
            _r = Residue(vm_object, name=_atom["resn"], index=_atom["resi"], chain=_chain)
            vm_object.residues[_atom["resi"]] = _r
            _chain.residues[_atom["resi"]] = _r
        _residue = _chain.residues[_atom["resi"]]
        
        atom = Atom(vm_object, name=_atom["name"], index=_atom["index"],
                    residue=_residue, chain=_chain, atom_id=atom_id,
                    occupancy=_atom["occupancy"], bfactor=_atom["bfactor"],
                    charge=_atom["charge"])
        atom.unique_id = unique_id
        atom._generate_atom_unique_color_id()
        _residue.atoms[atom_id] = atom
        vm_object.atoms[atom_id] = atom
        atom_id += 1
        unique_id += 1
    logger.debug("Time used to build the tree: {:>8.5f} secs".format(time.time() - initial))
    vm_object.frames = PDBFiles.get_coords_from_raw_frames(rawframes, atom_id, vismol_session.vm_config.n_proc)
    vm_object.mass_center = np.mean(vm_object.frames[0], axis=0)
    vm_object.find_bonded_and_nonbonded_atoms()
    return vm_object

def _load_psf_file(vismol_session, infile):
    """ Function doc
    """
    vismol_object = PSFFiles.load_PSF_topology_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_amber_top_file(vismol_session, infile):
    """ Function doc
    """
    vismol_object = AMBERFiles.load_amber_topology_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_xyz_file(vismol_session, infile):
    """ Function doc
    """
    vismol_object = XYZFiles.load_xyz_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def parse_file(vismol_session, infile):
    """ Function doc
    """
    show_molecule = True
    if infile[-3:] == "aux":
        vismol_object = _load_aux_file(vismol_session, infile)
    
    if infile[-3:] == "gro":
        vismol_object = _load_gro_file(vismol_session, infile)
    
    if infile[-4:] == "mol2":
        vismol_object = _load_mol2_file(vismol_session, infile)
    
    if infile[-3:] == "pdb":
        vismol_object = _load_pdb_file(vismol_session, infile)
    
    if infile[-3:] == "psf":
        vismol_object = _load_psf_file(vismol_session, infile)
        show_molecule = False
    
    if infile[-3:] == "top" or infile[-6:] == "prmtop":
        vismol_object = _load_amber_top_file(vismol_session, infile)
        show_molecule = False
    
    if infile[-3:] == "xyz":
        vismol_object = _load_xyz_file(vismol_session, infile)
    
    return vismol_object, show_molecule
