#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_object.py
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

import time
import numpy as np
from logging import getLogger
from vismol.model.atom import Atom
from vismol.model.bond import Bond
from vismol.model.chain import Chain
from vismol.model.residue import Residue
from vismol.model.molecular_properties import COLOR_PALETTE
from vismol.libgl.vismol_font import VismolFont
from vismol.libgl.representations import DotsRepresentation
from vismol.libgl.representations import LinesRepresentation
from vismol.libgl.representations import NonBondedRepresentation
from vismol.libgl.representations import PickingDotsRepresentation
from vismol.libgl.representations import ImpostorRepresentation
from vismol.libgl.representations import SticksRepresentation
from vismol.libgl.representations import SpheresRepresentation
from vismol.libgl.representations import DashedLinesRepresentation
# from vismol.libgl.representations import WiresRepresentation
# from vismol.libgl.representations import RibbonsRepresentation
import vismol.utils.c_distances as cdist

logger = getLogger(__name__)


class VismolObject:
    """ Visual Object contains the information necessary for openGL to draw 
        a model on the screen. Everything that is represented in the graphical 
        form is stored in the form of a VismolObject.
        
        Arguments
        
        name       = string  - Label that describes the object  
        atoms      = list of atoms  - [index, at_name, cov_rad,  at_pos, at_res_i, at_res_n, at_ch]
        vismol_session  = Vismol Session - Necessary to build the "atomtree_structure"
                     vismol_session contains the atom_id_counter (self.vm_session.atom_id_counter)
        trajectory = A list of coordinates - eg [ [x1,y1,z1, x2,y2,z2...], [x1,y1,z1, x2,y2,z2...]...]
                     One frame is is required at last.
    """
    def __init__(self, vismol_session, index, name="UNK", active=False, trajectory=None,
                 color_palette=None, bonds_pair_of_indexes=None):
        """ Class initialiser """
        self.vm_session = vismol_session
        self.vm_config = vismol_session.vm_config
        self.index = index         # import to find vboject in self.vm_session.vismol_objects_dic
        self.name = name
        self.active = active       # for "show and hide"   enable/disable
        self.frames = trajectory
        if color_palette is None:
            self.color_palette = COLOR_PALETTE[0]
        else:
            self.color_palette = color_palette #this is an integer num to access the color pick for carbon atoms (0 = green, 1 = purple, ...) 
        
        self.editing = False       # for translate and rotate  xyz coords 
        self.mass_center = np.zeros(3, dtype=np.float32)
        self.vm_font = VismolFont()
        self.coords = None
        
        self.atoms = {}
        self.residues = {}
        self.chains = {}
        self.atom_unique_id_dic = {}
        self.selected_atom_ids = set()
        self.bonds = None
        self.index_bonds = None # Pair of atoms, something like: [1, 3, 1, 17, 3, 4, 4, 20]
        self.non_bonded_atoms = None # Array of indexes
        self.cov_radii_array = None    # a list of covalent radius values for all  --> will be used to calculate de bonds
        
        self.representations = {}
        for rep_type in self.vm_config.representations_available:
            self.representations[rep_type] = None
        self.model_mat = np.identity(4, dtype=np.float32)
        self.trans_mat = np.identity(4, dtype=np.float32)
        
        self.core_representations = {"picking_dots":None, "picking_text":None}
        self.selection_dots_vao = None
        self.selection_dot_buffers = None
        self.picking_dots_vao = None
        self.picking_dot_buffers = None
        
        # self.dynamic_bons       = [] # Pair of atoms, something like: [0,1,1,2,3,4] 
        self.c_alpha_bonds = []
        self.c_alpha_atoms = []
    
    def build_core_representations(self):
        """ Function doc """
        self.core_representations["picking_dots"] = PickingDotsRepresentation(self,
                                                    self.vm_session.vm_glcore, active=True,
                                                    indexes=list(self.atoms.keys()))
        self.core_representations["dash"] = DashedLinesRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=self.index_bonds)
    
    def create_representation(self, rep_type="lines", indexes=None):
        """ Function doc """
        if rep_type == "dots":
            self.representations["dots"] = DotsRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=list(self.atoms.keys()))
        elif rep_type == "lines":
            self.representations["lines"] = LinesRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=self.index_bonds)
        elif rep_type == "nonbonded":
            self.representations["nonbonded"] = NonBondedRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=self.non_bonded_atoms)
        elif rep_type == "impostor":
            self.representations["impostor"] = ImpostorRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=list(self.atoms.keys()))
        elif rep_type == "sticks":
            self.representations["sticks"] = SticksRepresentation(self, self.vm_session.vm_glcore,
                                                                  active=True, indexes=self.index_bonds)
        elif rep_type == "spheres":
            self.representations["spheres"] = SpheresRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=list(self.atoms.keys()))
        
        # elif rep_type == "ribbons":
        #     self.representations["ribbons"] = RibbonsRepresentation(self, self.vm_session.vm_glcore,
        #                                                             active=True, indexes=indexes)
        # elif rep_type == "dotted_lines":
        #     self.representations["dotted_lines"] = LinesRepresentation(self, self.vm_session.vm_glcore,
        #                                                                active=True, indexes=indexes)
        else:
            logger.error("Representation {} not implemented".format(rep_type))
            raise NotImplementedError("Representation {} not implemented".format(rep_type))
    
    def _generate_color_vectors(self, colors_id_start, do_colors_raindow=True,
                                do_vdw_dot_sizes=True, do_cov_dot_sizes=True):
        """ (1) This method assigns to each atom of the system a 
            unique identifier based on the RGB color standard. 
            This identifier will be used in the selection function. 
            There are no two atoms with the same color ID in  
            
            (2) This method builds the "colors" np array that will 
            be sent to the GPU and which contains the RGB values 
            for each atom of the system.
        """
        atom_qtty = len(self.atoms)
        half = int(atom_qtty/2)
        quarter = int(atom_qtty/4)
        color_step = 1.0/(atom_qtty/4.0)
        red = 0.0
        green = 0.0
        blue = 1.0
        
        self.colors = np.empty([len(self.atoms), 3], dtype=np.float32)
        self.color_indexes = np.empty([len(self.atoms), 3], dtype=np.float32)
        if do_colors_raindow:
            self.color_rainbow = np.empty([len(self.atoms), 3], dtype=np.float32)
        if do_vdw_dot_sizes:
            self.vdw_dot_sizes = np.empty(len(self.atoms), dtype=np.float32)
        if do_cov_dot_sizes:
            self.cov_dot_sizes = np.empty(len(self.atoms), dtype=np.float32)
        
        for i, atom in self.atoms.items():
            self.colors[i] = atom.color
            self.color_indexes[i] = atom.color_id
            if do_vdw_dot_sizes: 
                self.vdw_dot_sizes[i] = atom.vdw_rad * 3
            if do_cov_dot_sizes: 
                self.cov_dot_sizes[i] = atom.cov_rad
            if do_colors_raindow:
                if i <= 1*quarter:
                    self.color_rainbow[i,:] = red, green, blue
                    green += color_step
                
                if (i >= 1*quarter) and (i <= 2*quarter):
                    self.color_rainbow[i,:] = red, green, blue
                    blue -= color_step
                
                if (i >= 2*quarter) and (i <= 3*quarter):
                    self.color_rainbow[i,:] = red, green, blue
                    red += color_step
                
                if (i >= 3*quarter) and (i <= 4*quarter):
                    self.color_rainbow[i,:] = red, green, blue
                    green -= color_step
    
    def find_bonded_and_nonbonded_atoms(self, selection=None, frame=0, gridsize=0.8,
                                         maxbond=2.4, tolerance=1.4):
        """ Function doc """
        if self.index_bonds is not None:
            logger.critical("It seems that there is already information about "\
                "the contacts in this VismolObject, trying to override the data "\
                "can produce serious problems :(")
        initial = time.time()
        atoms_frame_mask = np.zeros(len(self.atoms), dtype=np.bool)
        if self.cov_radii_array is None:
            self.cov_radii_array = np.empty(len(self.atoms), dtype=np.float32)
            for i, atom in self.atoms.items():
                self.cov_radii_array[i] = atom.cov_rad
        
        if selection is None:
            selection = self.atoms
            atoms_frame_mask[:] = True
        else:
            atoms_frame_mask[:] = False
            for atom in selection:
                atoms_frame_mask[atom.atom_id] = True
        
        cov_rads = self.cov_radii_array[atoms_frame_mask]
        coords = self.frames[frame][atoms_frame_mask]
        indexes = []
        gridpos_list = []
        for atom in selection.values():
            indexes.append(atom.atom_id)
            gridpos_list.append(atom.get_grid_position(gridsize=gridsize, frame=frame))
        logger.debug("Time used for preparing the atom mask, covalent radii list "\
                     "and grid positions: {}".format(time.time() - initial))
        
        initial = time.time()
        self.index_bonds = cdist.get_atomic_bonds_from_grid(indexes, coords,
                                        cov_rads, gridpos_list, gridsize, maxbond)
        msg = """Building grid elements  :
    Total number of Atoms   : {}
    Gridsize                : {}
    Bonds                   : {}
    Bonds calcultation time : {} seconds""".format(len(selection), gridsize,
                                len(self.index_bonds), time.time() - initial)
        logger.info(msg)
        self._bonds_from_pair_of_indexes_list()
        self._get_non_bonded_from_bonded_list()
    
    def _bonds_from_pair_of_indexes_list(self):
        """ Function doc 
        bonds_list = [[0,1] , [0,4] , [1,3], ...]
        """
        assert self.bonds is None
        self.bonds = []
        for i in range(0, len(self.index_bonds)-1, 2):
            index_i = self.index_bonds[i]
            index_j = self.index_bonds[i+1]
            bond = Bond(atom_i=self.atoms[index_i], atom_index_i=index_i,
                        atom_j=self.atoms[index_j], atom_index_j=index_j)
            self.bonds.append(bond)
            self.atoms[index_i].bonds.append(bond)
            self.atoms[index_j].bonds.append(bond)
        self.index_bonds = np.array(self.index_bonds, dtype=np.uint32)
    
    def _get_non_bonded_from_bonded_list(self):
        """ Function doc """
        assert self.non_bonded_atoms is None
        bonded_set = set(self.index_bonds)
        self.non_bonded_atoms = []
        for i, atom in self.atoms.items():
            if i in bonded_set:
                atom.nonbonded = False
            else:
                atom.nonbonded = True
                self.non_bonded_atoms.append(i)
        self.non_bonded_atoms.sort()
        self.non_bonded_atoms = np.array(self.non_bonded_atoms, dtype=np.int32)
    
    def _get_center_of_mass(self, frame=0):
        """ Function doc """
        if frame >= self.frames.shape[0]:
            logger.info("Frame {} is out of range for trajectory of size {}. \
                Using the last frame.".format(frame, self.frames.shape[0] - 1))
            frame = self.frames.shape[0] - 1
        return np.mean(self.frames[frame], axis=0)
    
    def set_model_matrix(self, mat):
        """ Function doc """
        self.model_mat = np.copy(mat)
    
    # def get_backbone_indexes(self):
    #     """ Function doc """
    #     chains_list = []
    #     bonds_pairs = []
    #     bonds_indexes = []
        
    #     self.c_alpha_bonds = []
    #     self.c_alpha_atoms = []
    #     for chain in self.chains:
    #         for residue in self.chains[chain].residues:
    #             if residue.is_protein:
    #                 for atom in residue.atoms:
    #                     if atom.name == "CA":
    #                         self.c_alpha_atoms.append(atom)
        
    #     for i in range(1, len(self.c_alpha_atoms)):
    #         atom_before = self.c_alpha_atoms[i-1]
    #         resi_before = atom_before.resi
    #         index_before = self.atoms.index(atom_before)
    #         atom = self.c_alpha_atoms[i]
    #         resi = atom.resi
    #         index = self.atoms.index(atom)
    #         if resi == resi_before + 1:
    #             bond = Bond(atom_i=atom_before, atom_index_i=index_before,
    #                         atom_j=atom, atom_index_j=index)
    #             distance = bond.distance()
    #             if distance < 4.0:
    #                 self.c_alpha_bonds.append(bond)
