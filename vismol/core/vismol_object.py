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
from model.atom import Atom
from model.bond import Bond
from model.chain import Chain
from model.residue import Residue
from model.molecular_properties import COLOR_PALETTE
from libgl.vismol_font import VismolFont
from libgl.representations import DotsRepresentation
from libgl.representations import LinesRepresentation
from libgl.representations import NonBondedRepresentation
# from libgl.representations import SticksRepresentation
# from libgl.representations import SpheresRepresentation
# from libgl.representations import ImpostorRepresentation
# from libgl.representations import WiresRepresentation
# from libgl.representations import RibbonsRepresentation
import utils.c_distances as cdist

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
        self.bonds = None
        self.index_bonds = None # Pair of atoms, something like: [1, 3, 1, 17, 3, 4, 4, 20]
        self.non_bonded_atoms = None # Array of indexes
        self.atoms_frame_mask = None
        self.cov_radii_array = None    # a list of covalent radius values for all  --> will be used to calculate de bonds
        
        self.representations = {}
        for rep_type in self.vm_config.representations_available:
            self.representations[rep_type] = None
        self.model_mat = np.identity(4, dtype=np.float32)
        self.trans_mat = np.identity(4, dtype=np.float32)
        
        self.atoms_by_chains    = {}    # A dictionary, access key is the chain id that connects with list of atom objects
        self.dynamic_bons       = [] # Pair of atoms, something like: [0,1,1,2,3,4] 
        self.noH_atoms = []
        self.c_alpha_bonds = []
        self.c_alpha_atoms = []
        self.residues_in_protein = []
        self.residues_in_solvent = []
        self.residues_ligands    = []
        self.atoms_in_protein    = [] # a list of atoms belonging to a protein
        self.atoms_in_solvent    = []
        self.selection_dots_vao    = None
        self.selection_dot_buffers = None
        self.target    = None
        self.unit_vec  = None
        self.distance  = None
        self.step      = None
        self.picking_dots_vao    = None
        self.picking_dot_buffers = None
        # """
        # This step is performed when no information about connections
        # between atoms is provided.
        # """
        # if auto_find_bonded_and_nonbonded:
        #     # this used just when the vobject is initialized
        #     self._find_bonded_and_nonbonded_atoms(selection=self.atoms, frame=0,
        #             gridsize=self.vm_session.vm_config.gl_parameters["gridsize"], 
        #             maxbond=self.vm_session.vm_config.gl_parameters["maxbond" ],
        #             tolerance=self.vm_session.vm_config.gl_parameters["bond_tolerance"])
        #     """ The nonbonded attribute of the atom object concerns representation.
        #         When true, I mean that the atom must be drawn with a small cross.
        #         You must assign the nonbonded attribute = True to atoms that are not bonded.
        #     """
        #     for index in self.non_bonded_atoms:
        #         self.atoms[index].nonbonded = True
        #     self._get_center_of_mass()
        
        # if bonds_pair_of_indexes:
        #     self._bonds_from_pair_of_indexes_list(bonds_pair_of_indexes)
        #     if len(self.non_bonded_atoms) == 0:
        #         self.import_non_bonded_atoms_from_bond()
    
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
        # elif rep_type == "sticks":
        #     self.representations["sticks"] = SticksRepresentation(self, self.vm_session.vm_glcore,
        #                                                           active=True, indexes=indexes)
        # elif rep_type == "ribbons":
        #     self.representations["ribbons"] = RibbonsRepresentation(self, self.vm_session.vm_glcore,
        #                                                             active=True, indexes=indexes)
        # elif rep_type == "spheres":
        #     self.representations["spheres"] = SpheresRepresentation(self, self.vm_session.vm_glcore,
        #                                                             active=True, indexes=indexes)
        #     self.representations["spheres"]._create_sphere_data()
        # elif rep_type == "dotted_lines":
        #     self.representations["dotted_lines"] = LinesRepresentation(self, self.vm_session.vm_glcore,
        #                                                                active=True, indexes=indexes)
        else:
            raise NotImplementedError("Representation {} not implemented".format(rep_type))
    
    def _generate_color_vectors(self, colors_id_start, do_colors=True,
                                do_colors_raindow=True, do_vdw_dot_sizes=True,
                                do_cov_dot_sizes=True):
        """ (1) This method assigns to each atom of the system a 
            unique identifier based on the RGB color standard. 
            This identifier will be used in the selection function. 
            There are no two atoms with the same color ID in  
            
            (2) This method builds the "colors" np array that will 
            be sent to the GPU and which contains the RGB values 
            for each atom of the system.
        """
        self.color_indexes = np.empty([len(self.atoms), 3], dtype=np.float32)
        self.colors = np.empty([len(self.atoms), 3], dtype=np.float32)
        for i, atom in self.atoms.items():
            self.colors[i] = atom.color
            self.color_indexes[i] = atom.color_id
        
        size = len(self.atoms)
        half = int(size/2)
        quarter = int(size/4)
        color_step = 1.0/(size/4)
        red = 0.0
        green = 0.0
        blue = 1.0
        
        # if do_colors:
        #     self.colors = []
        if do_colors_raindow:
            self.color_rainbow = []
        if do_vdw_dot_sizes:
            self.vdw_dot_sizes = []
        if do_cov_dot_sizes:
            self.cov_dot_sizes = []
        counter = 0
        for atom in self.atoms.values():
            #-------------------------------------------------------
            # (1)                  ID Colors
            #-------------------------------------------------------
            """
            i = atom.atom_id
            r = (i & 0x000000FF) >>  0
            g = (i & 0x0000FF00) >>  8
            b = (i & 0x00FF0000) >> 16
            """
            # if do_colors_idx:
            #     self.color_indexes.append(atom.color_id[0])
            #     self.color_indexes.append(atom.color_id[1])
            #     self.color_indexes.append(atom.color_id[2])
            """
            pickedID = r + g * 256 + b * 256*256
            atom.color_id = [r/255.0, g/255.0, b/255.0]
            #print (pickedID)
            self.vm_session.atom_dic_id[pickedID] = atom
            """
            #-------------------------------------------------------
            # (2)                   Colors
            #-------------------------------------------------------
            # if do_colors:
            #     self.colors.append(atom.color[0])
            #     self.colors.append(atom.color[1])
            #     self.colors.append(atom.color[2])
            #-------------------------------------------------------
            # (3)                  VdW list / cov_dot_sizes:
            #-------------------------------------------------------
            if do_vdw_dot_sizes: 
                self.vdw_dot_sizes.append(atom.vdw_rad * 3)
            if do_cov_dot_sizes: 
                self.cov_dot_sizes.append(atom.cov_rad)
            #-------------------------------------------------------
            # (4)                Rainbow colors
            #-------------------------------------------------------
            if do_colors_raindow:
                if counter <= 1*quarter:
                    self.color_rainbow.append(red)
                    self.color_rainbow.append(green)
                    self.color_rainbow.append(blue)
                    green += color_step
                
                if (counter >= 1*quarter) and (counter <= 2*quarter):
                    self.color_rainbow.append(red)
                    self.color_rainbow.append(green)
                    self.color_rainbow.append(blue)
                    blue -= color_step
                
                if (counter >= 2*quarter) and (counter <= 3*quarter):
                    self.color_rainbow.append(red)
                    self.color_rainbow.append(green)
                    self.color_rainbow.append(blue)
                    red += color_step
                
                if (counter >= 3*quarter) and (counter <= 4*quarter):
                    self.color_rainbow.append(red)
                    self.color_rainbow.append(green)
                    self.color_rainbow.append(blue)
                    green -= color_step
            counter += 1
        
        # if do_colors:
        #     self.colors = np.array(self.colors, dtype=np.float32)
        # if do_colors_idx:
        #     self.color_indexes = np.array(self.color_indexes, dtype=np.float32)
        if do_colors_raindow:
            self.colors_rainbow = np.array(self.color_rainbow, dtype=np.float32)
        if do_vdw_dot_sizes:
            self.vdw_dot_sizes = np.array(self.vdw_dot_sizes, dtype=np.float32)
        if do_cov_dot_sizes:
            self.cov_dot_sizes = np.array(self.cov_dot_sizes, dtype=np.float32)
        # print("-----------------Here-----------------", self.colors)
    
    def find_bonded_and_nonbonded_atoms(self, selection=None, frame=0, gridsize=1.33,
                                         maxbond=2.66, tolerance=1.4):
        """ Function doc """
        if self.index_bonds is not None:
            logger.critical("It seems that there is already information about \
                the contacts in this VismolObject, trying to override the data \
                can produce serious problems :(")
        initial = time.time()
        if self.atoms_frame_mask is None:
            self.atoms_frame_mask = np.zeros(len(self.atoms), dtype=np.bool)
        else:
            self.atoms_frame_mask[:] = False
        if self.cov_radii_array is None:
            self.cov_radii_array = np.empty(len(self.atoms), dtype=np.float32)
            for i, atom in self.atoms.items():
                self.cov_radii_array[i] = atom.cov_rad
        
        if selection is None:
            selection = self.atoms
            self.atoms_frame_mask[:] = True
        else:
            self.atoms_frame_mask[:] = False
            for atom in selection:
                self.atoms_frame_mask[atom.atom_id] = True
        
        cov_rads = self.cov_radii_array[self.atoms_frame_mask]
        coords = self.frames[frame][self.atoms_frame_mask]
        indexes = []
        gridpos_list = []
        for atom in selection.values():
            indexes.append(atom.atom_id)
            gridpos_list.append(atom.get_grid_position(gridsize=gridsize, frame=frame))
        logger.debug("Time used for preparing the atom mask, covalent radii list \
                     and grid positions: {}".format(time.time() - initial))
        
        initial = time.time()
        self.index_bonds = cdist.get_atomic_bonds_from_grid(indexes, coords,
                                        cov_rads, gridpos_list, gridsize, maxbond)
        logger.info("building grid elements  : {} secs".format(time.time() - initial))
        logger.info("Total number of Atoms   : {}".format(len(selection)))
        logger.info("Gridsize                : {}".format(gridsize))
        logger.info("Bonds                   : {}".format(len(self.index_bonds)))
        logger.info("Bonds calcultation time : {} secs".format(time.time() - initial))
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
            bond = Bond(atom_i=self.atoms[index_i],
                        atom_index_i=index_i,
                        atom_j=self.atoms[index_j],
                        atom_index_j=index_j)
            self.bonds.append(bond)
            # self.index_bonds.append(index_i)
            # self.index_bonds.append(index_j)
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
    
    def _generate_atomtree_structure(self, get_backbone_indexes=False):
        """ Function doc """
        print ("\n_generate_atomtree_structure starting")
        initial = time.time()
        self.atoms = []
        for raw_atom in self.raw_atoms_data:
            atom = Atom(index         = raw_atom["index"],
                        name          = raw_atom["name"],
                        symbol        = raw_atom["symbol"],
                        resi          = raw_atom["resi"],
                        resn          = raw_atom["resn"],
                        chain         = raw_atom["chain"],
                        atom_id       = self.vm_session.atom_id_counter,
                        occupancy     = raw_atom["occupancy"],
                        bfactor       = raw_atom["bfactor"],
                        charge        = raw_atom["charge"],
                        vismol_object = self)
            self.vm_session.atom_dic_id[self.vm_session.atom_id_counter] = atom
            self._add_new_atom_to_vobj(atom)
        self._get_center_of_mass()
        final = time.time()
        print ("_generate_atomtree_structure end -  total time: ", final - initial, "\n")
        if get_backbone_indexes:
            self.get_backbone_indexes()
        for chain in self.chains.keys():
            self.residues += self.chains[chain].residues
        return True
    
    def _add_new_atom_to_vobj(self, atom):
        """ Function doc """
        if atom.symbol == "H":
            pass
        else:
            self.noH_atoms.append(atom)
        
        if atom.chain in self.atoms_by_chains.keys():
            self.atoms_by_chains[atom.chain].append(atom)
        else:
            self.atoms_by_chains[atom.chain] = []
            self.atoms_by_chains[atom.chain].append(atom)
        
        if atom.chain in self.chains.keys():
            ch = self.chains[atom.chain]
        else:
            ch = Chain(name=atom.chain, label="UNK")
            self.chains[atom.chain] = ch
        """ This step checks if a residue has already been created and adds it
            to the respective chain.
        """
        if len(ch.residues) == 0:
            residue = Residue(name=atom.resn, index=atom.resi, chain=atom.chain,
                              vismol_object=self)
            atom.residue = residue
            residue.atoms.append(atom)
            ch.residues.append(residue)
            ch.residues_by_index[atom.resi] = residue
        elif atom.resi == ch.residues[-1].resi:
            atom.residue = ch.residues[-1]
            ch.residues[-1].atoms.append(atom)
        else:
            residue = Residue(name=atom.resn, index=atom.resi, chain=atom.chain,
                              vismol_object=self)
            atom.residue = residue
            residue.atoms.append(atom)
            ch.residues.append(residue)
            ch.residues_by_index[atom.resi] = residue
            #---------------------------------------------------------
            if residue.is_protein:
                self.residues_in_protein.append(residue)
            elif residue.is_solvent:
                self.residues_in_solvent.append(residue)
            else:
                self.residues_ligands.append(residue)
        
        if atom.name == "CA":
            ch.backbone.append(atom)
        self.atoms.append(atom)
        self.cov_radii_array.append(atom.cov_rad)
        self.vm_session.atom_dic_id[self.vm_session.atom_id_counter] = atom
        self.vm_session.atom_id_counter +=1
    
    def _get_center_of_mass(self, frame=0):
        """ Function doc
        """
        frame_size = len(self.frames) - 1
        
        if frame > frame_size:
            frame = frame_size
        
        # if len(self.noH_atoms) == 0:
        #     atoms = self.atoms
        # else:
        #     atoms = self.noH_atoms
        
        # com = np.zeros(3, dtype=np.float32)
        # initial  = time.time()
        # if len(self.frames) > 0:
        #     for atom in atoms:
        #         com += atom.coords(frame)
        # final = time.time()
        # self.mass_center = com / len(atoms)
        return np.mean(self.frames[frame], axis=0)
    
    def load_data_from_easyhybrid_serialization_file(self, d_atoms, frames, dynamic_bons):
        """ Function doc """
        print ("\nGenerate_chain_structure (easyhybrid_serialization_file) starting")
        initial           = time.time()
        self.frames       = frames
        self.atoms        = []
        self.dynamic_bons = dynamic_bons
        bonds_by_indexes = []
        for d_atom in d_atoms:
            for bond in d_atom["bonds"]:
                i = bond["atom_index_i"]
                j = bond["atom_index_j"]
                bonds_by_indexes.append([i,j])
            
            atom        = Atom(name          = d_atom["name"]                      ,
                               index         = d_atom["index"]                     ,
                               symbol        = d_atom["symbol"]                    , 
                               resi          = d_atom["resi"]                      ,
                               resn          = d_atom["resn"]                      ,
                               chain         = d_atom["chain"]                     ,
                               atom_id       = self.vm_session.atom_id_counter  , 
                               color         = d_atom["color"]                     , 
                               radius        = d_atom["radius"     ]               ,
                               vdw_rad       = d_atom["vdw_rad"    ]               ,
                               cov_rad       = d_atom["cov_rad"    ]               ,
                               ball_radius   = d_atom["ball_radius"]               ,
                               bonds_indexes = d_atom["bonds_indexes"]             ,
                               occupancy     = d_atom["occupancy"]                 ,
                               bfactor       = d_atom["bfactor"]                   ,
                               charge        = d_atom["charge"]                    ,
                               Vobject       = self                                ,
                               )
            atom.selected       = d_atom["selected"] 
            atom.lines          = d_atom["lines"] 
            atom.dots           = d_atom["dots"] 
            atom.nonbonded      = d_atom["nonbonded"] 
            atom.ribbons        = d_atom["ribbons"] 
            atom.ball_and_stick = d_atom["ball_and_stick"] 
            atom.sticks         = d_atom["sticks"] 
            atom.spheres        = d_atom["spheres"] 
            atom.surface        = d_atom["surface"] 
            atom.bonds_indexes  = d_atom["bonds_indexes"] 
            atom.bonds          = d_atom["bonds"] 
            atom.isfree         = d_atom["isfree"] 
            
            self.vm_session.atom_dic_id[self.vm_session.atom_id_counter] = atom
            self._add_new_atom_to_vobj(atom)  
        
        self._generate_color_vectors()
        self.bonds_from_pair_of_indexes_list(bonds_by_indexes)            
        if self.non_bonded_atoms == []:
            self.import_non_bonded_atoms_from_bond()
        self.get_backbone_indexes()
        final = time.time() 
        print ("_generate_atomtree_structure (easyhybrid_serialization_file) end -  total time: ", final - initial, "\n")
    
    def set_model_matrix(self, mat):
        """ Function doc
        """
        self.model_mat = np.copy(mat)
        return True
    
    def get_backbone_indexes(self):
        """ Function doc """
        chains_list = []
        bonds_pairs = []
        bonds_indexes = []
        
        self.c_alpha_bonds = []
        self.c_alpha_atoms = []
        for chain in self.chains:
            for residue in self.chains[chain].residues:
                if residue.is_protein:
                    for atom in residue.atoms:
                        if atom.name == "CA":
                            self.c_alpha_atoms.append(atom)
        
        for i in range(1, len(self.c_alpha_atoms)):
            atom_before = self.c_alpha_atoms[i-1]
            resi_before = atom_before.resi
            index_before = self.atoms.index(atom_before)
            atom = self.c_alpha_atoms[i]
            resi = atom.resi
            index = self.atoms.index(atom)
            if resi == resi_before + 1:
                bond = Bond(atom_i=atom_before, atom_index_i=index_before,
                            atom_j=atom, atom_index_j=index)
                distance = bond.distance()
                if distance < 4.0:
                    self.c_alpha_bonds.append(bond)
    
    def find_bonded_and_nonbonded_by_selection(self, selection=None, frame=0,
                                               gridsize=1.4, tolerance=1.4, maxbond=2.6):
        """ Function doc """
        if selection == None:
            selection = self.atoms
        atoms_list = []
        indexes = []
        cov_rad = np.array(self.cov_radii_array, dtype=np.float32)
        coords = self.frames[frame]
        gridpos_list = []
        for atom in selection:
            indexes.append(atom.index-1)
            gridpos_list.append(atom.get_grid_position(gridsize=gridsize, frame=frame))
        bonds_pair_of_indexes = cdist.get_atomic_bonds_from_grid(indexes, coords,
                                        cov_rad, gridpos_list, gridsize, maxbond)
        return  bonds_pair_of_indexes 
