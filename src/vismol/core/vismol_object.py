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
from vismol.model.molecule import Molecule
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
from vismol.libgl.representations import LabelRepresentation
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
        self.molecules = {}

        self.atom_unique_id_dic = {}
        self.selected_atom_ids = set()
        self.bonds = None
        self.index_bonds = None # Pair of atoms, something like: [1, 3, 1, 17, 3, 4, 4, 20]
        self.non_bonded_atoms = None # Array of indexes
        self.cov_radii_array = None    # a list of covalent radius values for all  --> will be used to calculate de bonds
        
          
        self.topology = {} # {92: [93, 99], 93: [92, 94, 96, 100], 99: [92],...}
                           # important to define molecules
                           
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
        
        self.dynamic_bonds  = [] # Pair of atoms, something like: [[0,1,1,2,3,4] , [0,1,1,2], ...]
                                # Like self.index_bonds but for each frame
        self.c_alpha_bonds = []
        self.c_alpha_atoms = []
    
        ''' Cell and Symmetry'''
        self.cell_parameters  = None
        self.cell_coordinates = None
        self.cell_indexes     = None
        self.cell_colors      = None
        self.cell_bonds       = None
        self.representations['labels'] = None
    
    
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
                                                    active=True, indexes=list(self.atoms.keys()) )
        
        elif rep_type == "picking_spheres":
            self.representations["picking_spheres"] = SpheresRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=list(self.atoms.keys()), mode =1 )
        
        elif rep_type == "vdw_spheres":
            self.representations["vdw_spheres"] = SpheresRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=list(self.atoms.keys()), vdw =True)

        elif rep_type == "dash":
            self.representations["dash"] = DashedLinesRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=self.index_bonds)
        elif rep_type == "ribbons":
            self.representations["ribbons"] = SticksRepresentation(self, self.vm_session.vm_glcore,
                                                                    active=True, indexes=self.index_bonds, name  = 'ribbons')
        elif rep_type == "dynamic":
            #print(self.dynamic_bonds)
            self.representations["dynamic"] = SticksRepresentation(self, self.vm_session.vm_glcore,
                                                                  active=True, indexes=self.index_bonds, is_dynamic = True)
        
        elif rep_type == "labels":
            #print(self.dynamic_bonds)
            self.representations["labels"] = LabelRepresentation(vismol_object  = self  ,  
                                                                  vismol_glcore = self.vm_session.vm_glcore , 
                                                                  indexes       = [0,1,2] , 
                                                                  labels        = None     , 
                                                                  color         = [1, 1, 0, 1])
        
        
        
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
    
    def define_bonds_from_external (self, index_bonds = [], internal = True):
        """ Function doc """
        if internal is True:
            self.index_bonds = index_bonds

            self._bonds_from_pair_of_indexes_list()
            self._get_non_bonded_from_bonded_list()
            
            self._generate_topology_from_index_bonds()
            self.define_molecules()
            self.define_Calpha_backbone()
        else:
            return index_bonds
    
    def find_bonded_and_nonbonded_atoms(self, selection=None, frame=0, gridsize=0.8,
                                         maxbond=2.4, tolerance=1.4, internal = True):
        """
            Receives a dictionary as a selection:

            selection = {'atom_id': atom_object, ...} ()

            frame - (frame number)

            grid_size - (usually should not be changed. Default is: 0.8)

            maxbond - (Size of the maximum bond to be monitored. Bonds greater than "maxbond" 
            may be disregarded. Default is: 2.4 A)

            tolerance (Safety factor that multiplies the term (ra_cov + rb_cov)**2. Default is: 1.4)

            internal (If True, calculates the bindings for the object itself. Used in the object's 
            genesis. If False, returns a list of atom_ids like  [0,1, # bond  between 0 and 1 
                                                                 1,2, # bond  between 1 and 2
                                                                 0,3] # bond  between 0 and 3 )
        """
        if internal:
            if self.index_bonds is not None:
                logger.critical("It seems that there is already information about "\
                    "the contacts in this VismolObject, trying to override the data "\
                    "can produce serious problems :(")
        
        
        initial = time.time()
        atoms_frame_mask = np.zeros(len(self.atoms), bool)
        if self.cov_radii_array is None:
            self.cov_radii_array = np.empty(len(self.atoms), dtype=np.float32)
            for i, atom in self.atoms.items():
                self.cov_radii_array[i] = atom.cov_rad
        
        if selection is None:
            selection = self.atoms
            atoms_frame_mask[:] = True
        else:
            atoms_frame_mask[:] = False
            for atom in selection.values():
                atoms_frame_mask[atom.atom_id] = True
        
        #cov_rads = self.cov_radii_array[atoms_frame_mask]
        #coords = self.frames[frame][atoms_frame_mask]
        cov_rads = self.cov_radii_array
        coords = self.frames[frame]#[atoms_frame_mask]
        
        indexes = []
        gridpos_list = []
        for atom in selection.values():
            indexes.append(atom.atom_id)
            gridpos_list.append(atom.get_grid_position(gridsize=gridsize, frame=frame))
        logger.debug("Time used for preparing the atom mask, covalent radii list "\
                     "and grid positions: {}".format(time.time() - initial))
        if internal is True:
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
            
            #bachega's code
            initial = time.time()
            self._generate_topology_from_index_bonds()
            
            
            #try:
            self.define_molecules()
            #except:
            #    print('Error: failure when trying to find molecules')
            
            
            final = time.time()
            print('        Defining molecule indexes: ', final - initial)
            #print ('\n',self.topology, '\n')
            self.define_Calpha_backbone()
        else:
            index_bonds = cdist.get_atomic_bonds_from_grid(indexes, coords,
                                            cov_rads, gridpos_list, gridsize, maxbond)
            return index_bonds
    
    def _bonds_from_pair_of_indexes_list(self):
        """ Function doc 
        self.index_bonds = [0,1  , 0,4  ,  1,3  , ...]
        self.bonds = [bond1(obj), bond2(obj), bond3(obj), ...] 
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
    
    def _generate_topology_from_index_bonds(self, bonds = None):
        """ 
        bonds = [92,93  ,  92,99  ,  ...]
        
        Returns a graph in dictionary form (this may in turn be 
        needed to determine which objects are molecules).
        
        defines: self.topology = {92: [93, 99], 93: [92, 94, 96, 100], 99: [92],...}
        
        """
        if bonds is None:
            bonds =  self.index_bonds
        else:
            pass

        
        bonds_pairs = []
        topology    = {}
        
        for i  in range(0, len(bonds),2):
            bonds_pairs.append([bonds[i], bonds[i+1]])
        
        
        for bond in bonds_pairs:
            if bond[0] in topology.keys():
                topology[bond[0]].append(bond[1])
            else:
                topology[bond[0]] = []
                topology[bond[0]].append(bond[1])
            
            
            if bond[1] in topology.keys():
                topology[bond[1]].append(bond[0])
            else:
                topology[bond[1]] = []
                topology[bond[1]].append(bond[0])
        self.topology = topology

    def define_molecules (self):
        """ Function doc 
        self.topology  = It is a graph written in the form of a dictionary:
                        {
                         index1 : [index2 , index3, ...], 
                         index2 : [index1 , index4, ...]
                         ...}
        
        groups = [{0, 1, 2}, {3, 4, 5, 6}, {8, 9, 7}]
        
        self.atoms =  {
                       0: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323cbcdf0>, 
                       1: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323cbcfd0>, 
                       2: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323cbcf40>, 
                       3: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323e38250>, 
                       4: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323e38040>, 
                       5: ...
                       }
        
        This function populates the molecule dictionary (self.molecules), each 
        molecule object is created very similarly to the residue object
        
        """
        #try:
        #groups = find_groups(self.topology)
        groups = find_connected_components(self.topology)
        mol_index = 0
        for mol_index, group in enumerate( groups ):
            atoms = {}
            molecule = Molecule(self, name="UNK", index = mol_index)
            
            for atom_index in list(group):
                molecule.atoms[atom_index] = self.atoms[atom_index]
                self.atoms[atom_index].molecule = molecule
                
            self.molecules[mol_index] = molecule
        
        #-------------------------------------------------------------
        # non_bonded_atoms
        # Should be here, ohterwise we will have selection problems
        #-------------------------------------------------------------

        for atom_index in self.non_bonded_atoms:
            #try:
            mol_index += 1
            molecule = Molecule(self, name="UNK", index = mol_index)
            molecule.atoms[atom_index] = self.atoms[atom_index]
            self.atoms[atom_index].molecule = molecule
            self.molecules[mol_index] = molecule

            #print(self.molecules)
        #except:
        #    print('Failure to determine the list of molecules!')
            
    def define_Calpha_backbone (self):
        """ Function doc 
        Verifica quais conexões entre c_alphas são válidas.
        """
        
        
        self.c_alpha_bonds = []
        self.c_alpha_atoms = []
        
        #
        # Building the self.c_alpha_atoms dict
        # {atom_id : atom_object, ...}
        #
        for c_index, chain in self.chains.items():
            for r_index, residue in chain.residues.items():
                if residue.is_protein:
                    for a_index, atom in residue.atoms.items():
                        if atom.name == "CA":
                            self.c_alpha_atoms.append(atom)

        
        for i in range(1, len(self.c_alpha_atoms)):

            atom_before  = self.c_alpha_atoms[i-1]
            resi_before  = atom_before.residue.index
            index_before = atom_before.atom_id
            

            atom  = self.c_alpha_atoms[i]
            resi  = atom.residue.index
            index = atom.atom_id
            
            '''Checks whether the two residues are in sequence 
            (otherwise there is a break in the backbone structure)'''
            if resi == resi_before + 1:
                
                bond = Bond(atom_i=atom_before, atom_index_i=index_before,
                            atom_j=atom, atom_index_j=index)
                
                distance = bond.distance()
                if distance < 4.0:
                    self.c_alpha_bonds.append(bond)
        
        #print(self.c_alpha_bonds)
        #print(self.c_alpha_atoms)
        

    def _calculate_unit_cell_vertices(self, a, b, c, alpha, beta, gamma):
        '''
        Returns the list of vertices positions of 
        a box with parameters a, b, c, alpha, beta, gamma.
        '''
        
        # Convert angles to radians
        alpha_rad = np.radians(alpha)
        beta_rad  = np.radians(beta)
        gamma_rad = np.radians(gamma)

        # Calculate unit cell vectors
        v1 = np.array([a, 0, 0])
        v2 = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
        v3_x = c * np.cos(beta_rad)
        v3_y = c * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)
        v3_z = np.sqrt(c**2 - v3_x**2 - v3_y**2)
        v3 = np.array([v3_x, v3_y, v3_z])

        # Define the eight vertices of the unit cell
        vertices = [
            np.array([0, 0, 0]),
            v1,
            v2,
            v2 + v1,
            v3,
            v3 + v1,
            v3 + v2,
            v3 + v2 + v1
        ]

        return vertices
    
    def set_cell (self, a, b, c, alpha, beta, gamma, color = [0.5, 0.5, 0.5]):
        """
        Assign the cell parameters to the 
        "vismol_object" object.
        
        Arguments a, b, c, alpha, beta, gamma 
        must be obtained externally.
        
        """
        
        self.cell_parameters = {'a'     : a, 
                                'b'     : b, 
                                'c'     : c, 
                                'alpha' : alpha, 
                                'beta'  : beta, 
                                'gamma' : gamma
                               }
         
        vertices = self._calculate_unit_cell_vertices(a, b, c, alpha, beta, gamma)


        '''
        The coordinates follow the same structure as the atoms, 
        they are organized in a trajectory structure with only one 
        frame.
        '''
        self.cell_coordinates = np.empty([1, 8, 3], dtype=np.float32)
        self.cell_indexes     = []
       
        
        '''
        The connections between the vertices of the boxes 
        are pre-established.
        '''
        self.cell_bonds       = [0,1, 0,2, 2,3, 0,4 , 1,3 , 1,5 , 2,3 , 2,6 , 3,7 , 4,6 , 4,5, 5,7, 6,7]       
        self.cell_bonds       = np.array(self.cell_bonds, dtype=np.uint32)
        
        
        '''
        The colors of the vertices follow the same 
        structure used in atoms
        '''
        self.cell_colors      = np.empty([8, 3], dtype=np.float32)

        for i, vertex in enumerate(vertices, 0):
        
            print("Vertex {}: {:7.3f} {:7.3f} {:7.3f}".format(i, vertex[0] ,vertex[1] ,vertex[2]))
            self.cell_indexes.append(i)
            
            self.cell_colors[i] = np.array(color, dtype=np.float32)
            self.cell_coordinates[0,i,:] = vertex[0] ,vertex[1] ,vertex[2]  
            

        print('\n\n')
        print (self.cell_colors)







def find_connected_components(graph):
    """
    
    graph = topology (self.topology)
    
    eg: 

        self.topology = {0: [1, 2, 1], 1: [0, 0], 2: [0], 3: [6, 4, 5], 6: [3], 4: [3], 5: [3], 7: [8, 9], 8: [7], 9: [7]}
    
    
    
    This version of the function also takes a graph represented as a dictionary and 
    returns a list of connected components, where each connected component is a 
    list of nodes.

    The function starts by initializing an empty set called visited to keep track 
    of the nodes that have been visited and an empty list called components to 
    store the connected components.

    The function then iterates over each node in the graph and checks if it has 
    been visited yet. If the node has not been visited, the function starts a DFS 
    from that node.

    The DFS is implemented using a stack-based approach. The function initializes 
    a stack with the starting node and an empty list called component to store the 
    nodes in the connected component.

    The function then enters a loop that continues as long as the stack is not 
    empty. In each iteration of the loop, the function pops a node from the stack 
    and checks if it has been visited yet. If the node has not been visited, the 
    function adds it to the visited set, appends it to the component list, and 
    adds its unvisited neighbors to the stack.

    After the DFS has completed, the function appends the component list to the 
    components list.

    Finally, the function returns the components list, which contains the connected 
    components of the graph. Each connected component is represented as a list of 
    nodes.
    
    eg:
        components = [[0, 1, 2], [3, 5, 4, 6], [7, 9, 8]]
    
    """
    
    
    visited = set()
    components = []

    for start_node in graph:
        if start_node not in visited:
            stack = [start_node]
            component = []
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    component.append(node)
                    stack.extend([neighbor for neighbor in graph[node] if neighbor not in visited])
            components.append(component)
    return components

def DFS(graph, node, visited):
    ''' 
        The DFS function takes a graph, a node, and a set of 
        visited nodes as inputs, and performs a depth-first 
        search starting from the node. 
        
        is not being used

    '''
    visited.add(node)
    for neighbor in graph[node]:
        if neighbor not in visited:
            DFS(graph, neighbor, visited)

def find_groups(graph):
    '''
    
    is not being used due: 
    
    Python: maximum recursion depth exceeded while calling a Python object

    for more info access: 
    https://stackoverflow.com/questions/6809402/python-maximum-recursion-depth-exceeded-while-calling-a-python-object
        
    '''
    visited = set()
    groups = []
    for node in graph:
        if node not in visited:
            group = set()
            DFS(graph, node, group)
            groups.append(group)
            visited |= group
    return groups
