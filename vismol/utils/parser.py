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

from utils import AUXFiles
from utils import GROFiles
from utils import MOL2Files
from utils import PDBFiles
from utils import PSFFiles
from utils import AMBERFiles
from utils import XYZFiles


def _load_aux_file(infile, vismol_session):
    """ Function doc """
    vismol_object = AUXFiles.load_aux_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_gro_file(infile, vismol_session):
    """ Function doc
    """
    vismol_object = GROFiles.load_gro_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_mol2_file(infile, vismol_session):
    """ Function doc """
    vismol_object = MOL2Files.load_mol2_files(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_pdb_file(infile, vismol_session):
    """ Function doc
    """
    vismol_object = PDBFiles.load_pdb_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_psf_file(infile, vismol_session):
    """ Function doc
    """
    vismol_object = PSFFiles.load_PSF_topology_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_amber_top_file(infile, vismol_session):
    """ Function doc
    """
    vismol_object = AMBERFiles.load_amber_topology_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def _load_xyz_file(infile, vismol_session):
    """ Function doc
    """
    vismol_object = XYZFiles.load_xyz_file(infile, vismol_session)
    vismol_object.set_model_matrix(vismol_session.vm_widget.vm_glcore.model_mat)
    return vismol_object

def parse_file(infile, vismol_session):
    """ Function doc
    """
    show_molecule = True
    if infile[-3:] == "aux":
        vismol_object = _load_aux_file(infile, vismol_session)
    
    if infile[-3:] == "gro":
        vismol_object = _load_gro_file(infile, vismol_session)
    
    if infile[-4:] == "mol2":
        vismol_object = _load_mol2_file(infile, vismol_session)
    
    if infile[-3:] == "pdb":
        vismol_object = _load_pdb_file(infile, vismol_session)
    
    if infile[-3:] == "psf":
        vismol_object = _load_psf_file(infile, vismol_session)
        show_molecule = False
    
    if infile[-3:] == "top" or infile[-6:] == "prmtop":
        vismol_object = _load_amber_top_file(infile, vismol_session)
        show_molecule = False
    
    if infile[-3:] == "xyz":
        vismol_object = _load_xyz_file(infile, vismol_session)
    
    return vismol_object, show_molecule
