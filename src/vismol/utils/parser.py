#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import time
import numpy as np
from typing import Tuple
from logging import getLogger
from vismol.utils import AUXFiles
from vismol.utils import GROFiles
from vismol.utils import MOL2Files
from vismol.utils import PDBFiles
from vismol.utils import PSFFiles
from vismol.utils import AMBERFiles
from vismol.utils import XYZFiles


logger = getLogger(__name__)


def _load_aux_file(vm_session: "VismolSession", infile: str) -> "VismolObject":
    """
    Load an auxiliary file in the VismolSession.
    
    Args:
        vm_session (VismolSession): The VismolSession object.
        infile (str): The path to the auxiliary file.
    
    Returns:
        VismolObject: The loaded VismolObject.
    
    """
    vm_object = AUXFiles.load_aux_file(vm_session, infile)
    # vm_object.set_model_matrix(vm_session.vm_widget.vm_glcore.model_mat)
    return vm_object


def _load_gro_file(vm_session: "VismolSession", infile: str) -> "VismolObject":
    """
    Load a GRO file in the VismolSession.
    
    Args:
        vm_session (VismolSession): The VismolSession object.
        infile (str): The path to the GRO file.
    
    Returns:
        VismolObject: The loaded VismolObject.
    
    """
    vm_object = GROFiles.load_gro_file(vm_session, infile)
    # vm_object.set_model_matrix(vm_session.vm_widget.vm_glcore.model_mat)
    return vm_object


def _load_mol2_file(vm_session: "VismolSession", infile: str) -> "VismolObject":
    """
    Load a MOL2 file in the VismolSession.
    
    Args:
        vm_session (VismolSession): The VismolSession object.
        infile (str): The path to the MOL2 file.
    
    Returns:
        VismolObject: The loaded VismolObject.
    
    """
    vm_object = MOL2Files.load_mol2_files(vm_session, infile)
    # vm_object.set_model_matrix(vm_session.vm_widget.vm_glcore.model_mat)
    return vm_object


def _load_pdb_file(vm_session: "VismolSession", infile: str) -> "VismolObject":
    """
    Load a PDB file in the VismolSession.
    
    Args:
        vm_session (VismolSession): The VismolSession object.
        infile (str): The path to the PDB file.
    
    Returns:
        VismolObject: The loaded VismolObject.
    
    """
    vm_object = PDBFiles.load_pdb_file(vm_session, infile)
    # vm_object.set_model_matrix(vm_session.vm_widget.vm_glcore.model_mat)
    return vm_object


def _load_psf_file(vm_session: "VismolSession", infile: str) -> "VismolObject":
    """
    Load a PSF file in the VismolSession.
    
    Args:
        vm_session (VismolSession): The VismolSession object.
        infile (str): The path to the PSF file.
    
    Returns:
        VismolObject: The loaded VismolObject.
    
    """
    vm_object = PSFFiles.load_PSF_topology_file(vm_session, infile)
    # vm_object.set_model_matrix(vm_session.vm_widget.vm_glcore.model_mat)
    return vm_object


def _load_amber_top_file(vm_session: "VismolSession", infile: str) -> "VismolObject":
    """
    Load an AMBER topology file in the VismolSession.
    
    Args:
        vm_session (VismolSession): The VismolSession object.
        infile (str): The path to the AMBER topology file.
    
    Returns:
        VismolObject: The loaded VismolObject.
    
    """
    vm_object = AMBERFiles.load_amber_topology_file(vm_session, infile)
    # vm_object.set_model_matrix(vm_session.vm_widget.vm_glcore.model_mat)
    return vm_object


def _load_xyz_file(vm_session: "VismolSession", infile: str) -> "VismolObject":
    """
    Load an XYZ file in the VismolSession.
    
    Args:
        vm_session (VismolSession): The VismolSession object.
        infile (str): The path to the XYZ file.
    
    Returns:
        VismolObject: The loaded VismolObject.
    
    """
    vm_object = XYZFiles.load_xyz_file(vm_session, infile)
    # vm_object.set_model_matrix(vm_session.vm_widget.vm_glcore.model_mat)
    return vm_object


def parse_file(vm_session: "VismolSession", infile: str) -> Tuple["VismolObject", bool]:
    """
    Parse a file and load the corresponding VismolObject.
    
    Args:
        vm_session (VismolSession): The VismolSession object.
        infile (str): The path to the input file.
    
    Returns:
        Tuple[VismolObject, bool]: The loaded VismolObject and a boolean
                                   indicating whether to show the molecule or not.
    
    """
    show_molecule = True
    
    if infile.endswith("aux"):
        vm_object = _load_aux_file(vm_session, infile)
    elif infile.endswith("gro"):
        vm_object = _load_gro_file(vm_session, infile)
    elif infile.endswith("mol2"):
        vm_object = _load_mol2_file(vm_session, infile)
    elif infile.endswith("pdb"):
        vm_object = _load_pdb_file(vm_session, infile)
    elif infile.endswith("psf"):
        vm_object = _load_psf_file(vm_session, infile)
        show_molecule = False
    elif infile.endswith("top") or infile.endswith("parm7"):
        vm_object = _load_amber_top_file(vm_session, infile)
        show_molecule = False
    elif infile.endswith("xyz"):
        vm_object = _load_xyz_file(vm_session, infile)
    else:
        raise ValueError(f"Unsupported file format: {infile}")
    
    return vm_object, show_molecule

