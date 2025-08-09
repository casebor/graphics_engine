#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import time
import numpy as np
from logging import getLogger
from vismol.model.bond import Bond
import vismol.utils.c_distances as cdist


logger = getLogger(__name__)


class Molecule:
    """
    The Molecule contains all the molecular information, i.e., the topology, the
    atoms, bonds, etc.
    
    Attributes:
        vismol_object (VismolObject): the Vismol object that this molecule
                belongs to. The Vismol object is a connection between the
                molecular structure and the OpenGL engine.
        name (str): The name of the molecule. Not sure where or how is this
                attribute required.
        trajectory (np.array): An array of coordinates representing the
                trajectory frames. Each frame is represented as an array of
                [x, y, z] coordinates (optional). The final array shape is 
                [frames, atoms, 3]
        chains (Dict): Dictionary containing the chains of the molecule. A molecule
                can contain only chains, not residues, however, the atoms are
                present at this level as well because they are used for the
                selection in a fast manner. TODO: This can be improved.
        geom_center (np.array): A Numpy array containing geometry center of the
                molecule.
        atoms (Dict): Dictionary containing all the atoms of the molecule.
        atom_unique_id_dic (Dict): Dictionary containing the color coding for
                each atom in the VisMol object used for selection. This code id
                is unique. TO_REVIEW
        selected_atom_ids (Set): Set to store the IDs of the selected atoms. TO_REVIEW
        bonds (List): List of Bond objects used represent molecular bonds
                between atoms. TO_REVIEW
        index_bonds (List): List of indexes of atoms that make the bonds. TO_REVIEW
        non_bonded_atoms (List): List of indexes of atoms that do not make
                any bonds. TO_REVIEW
        cov_radii_array
    """
    
    def __init__(self, vismol_object: "VismolObject", name: str="UNK",
                 trajectory: np.array=None):
        """
        Args:
            vismol_object (VismolObject): the Vismol object that this molecule
                    belongs to. The Vismol object is a connection between the
                    molecular structure and the OpenGL engine.
            name (str): The name of the molecule. Not sure where or how is this
                    attribute required.
            trajectory (np.array): An array of coordinates representing the
                    trajectory frames. Each frame is represented as an array of
                    [x, y, z] coordinates (optional). The final array shape is 
                    [frames, atoms, 3]
        """
        self.vm_object = vismol_object
        self.name = name
        self.frames = trajectory
        self.chains = {}
        self.geom_center = np.zeros(3, dtype=np.float32)
        self.atoms = {}
        self.atom_unique_id_dic = {}
        
        
        self.bonds = None       # Bond objects representing connections between atoms
        self.index_bonds = None # Pair of atoms, something like: [1, 3, 1, 17, 3, 4, 4, 20]
                                # Pair of atoms used to define bonds (set as None if not provided)
        
        self.non_bonded_atoms = None # Array of indexes
                                     # Array of indexes for non-bonded atoms (not yet defined in the code)
                                     
        self.cov_radii_array = None  # a list of covalent radius values for all  --> will be used to calculate de bonds
                                     # List of covalent radius values for all atoms (not yet defined in the code)
        
          
        self.topology = {} # {92: [93, 99], 93: [92, 94, 96, 100], 99: [92],...}
                           # important to define molecules
                           
    
    def _is_protein(self):
        """ Function doc """
        # is it a protein residue?
        if self.name in residues_dictionary.keys():
            self.is_protein = True
        # is it a salvent molecule?
        if self.name in solvent_dictionary.keys():
            self.is_solvent = True
    
    def geometry_center(self, frame=0):
        """ Function doc """
        if frame > len(self.vm_object.frames)-1:
            frame = len(self.vm_object.frames)-1
        gc = np.zeros(3, dtype=np.float32)
        for atom in self.atoms.values():
            gc += atom.coords(frame)
        gc /= len(self.atoms.values())
        return gc
    
    def get_center_of_mass(self, mass=False, frame=0):
        """ Function doc """
        frame_size = len(self.vm_object.frames)-1
        
        if frame <= frame_size:
            pass
        else:
            frame = frame_size
        
        total = len(self.atoms)
        sum_x = 0.0
        sum_y = 0.0
        sum_z = 0.0
        
        for atom in self.atoms.values():
            coord = atom.coords (frame)
            sum_x += coord[0]
            sum_y += coord[1]
            sum_z += coord[2]
        
        self.geom_center = np.array([sum_x / total,
                                     sum_y / total, 
                                     sum_z / total])
    
    def _get_center_of_mass(self, frame=0):
        """ Function doc """
        if frame >= self.frames.shape[0]:
            logger.info("Frame {} is out of range for trajectory of size {}. \
                Using the last frame.".format(frame, self.frames.shape[0] - 1))
            frame = self.frames.shape[0] - 1
        return np.mean(self.frames[frame], axis=0)
    
    def define_bonds_from_external(self, index_bonds = [], internal = True):
        """ Function doc """
        if internal is True:
            self.index_bonds = index_bonds

            self._bonds_from_pair_of_indexes_list()
            self._get_non_bonded_from_bonded_list()
            
            self._generate_topology_from_index_bonds()
            # self.define_molecules()
            self.define_Calpha_backbone()
        else:
            return index_bonds
    
    def build_bonded_and_nonbonded_atoms(self, frame=0, gridsize=1.2,
                                        maxbond=2.4, tolerance=1.4):
        """
        Function doc: Determines bonded and nonbonded atoms based on selection in a VismolObject.
        
        Parameters:
        - selection: A dictionary containing the selected atoms (optional).
        - frame: Frame number for the atom coordinates (default is 0).
        - gridsize: Grid size for bond calculation (usually should not be changed, default is 0.8).
        - maxbond: Size of the maximum bond to be monitored (bonds greater than maxbond may be disregarded, default is 2.4 Å).
        - tolerance: Safety factor that multiplies the term (ra_cov + rb_cov)**2 (default is 1.4).
        - internal: If True, calculates the bindings for the object itself (used in the object's genesis).
                    If False, returns a list of atom_ids representing the bonds between atoms.
        
        Returns:
        If internal is True, the function updates the VismolObject's internal attributes:
        - index_bonds: List of pairs of atom indexes representing the bonded atoms.
        - bonds: List of Bond objects representing the bonds between atoms.
        - topology: Dictionary representing the atom topology and connectivity.
        - Non-bonded interactions and other attributes are also updated.
        
        If internal is False, the function returns a list of atom_ids representing the bonds between atoms.
        
        Note: This function calculates bonds between atoms based on their positions and covalent radii.
        
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
        # Check if internal is True and there is already information about contacts.
        if self.index_bonds is not None:
            logger.critical("It seems that there is already information about "\
                "the contacts in this VismolObject, trying to override the data "\
                "can produce serious problems :(")
        
        # Initialize variables and data structures
        initial = time.time()
        if self.cov_radii_array is None:
            self.cov_radii_array = np.empty(len(self.atoms), dtype=np.float32)
            for i, atom in self.atoms.items():
                self.cov_radii_array[i] = atom.cov_rad
        
        # Extract relevant data for bond calculation
        cov_rads = self.cov_radii_array
        coords = self.frames[frame]
        
        # Generate indexes and grid positions for each atom in the selection
        indexes = []
        gridpos_list = []
        for atom in self.atoms.values():
            indexes.append(atom.atom_id)
            gridpos_list.append(atom.get_grid_position(gridsize=gridsize, frame=frame))
        logger.debug("Time used for preparing the atom mask, covalent radii list "\
                         "and grid positions: {}".format(time.time() - initial))
        
        # Calculate bonds based on grid positions and covalent radii
        initial = time.time()
        self.index_bonds = cdist.get_atomic_bonds_from_grid(indexes, coords,
                                        cov_rads, gridpos_list, gridsize,
                                        maxbond, tolerance)
        msg = """Building grid elements  :
    Total number of Atoms   : {}
    Gridsize                : {}
    Bonds                   : {}
    Bonds calcultation time : {} seconds""".format(len(self.atoms), gridsize,
                                len(self.index_bonds), time.time() - initial)
        logger.info(msg)
        
        # Create Bond objects and update atom bond lists
        self._bonds_from_pair_of_indexes_list()
        
        # Generate non-bonded interactions from bonded atoms
        self._get_non_bonded_from_bonded_list()
        
        # Generate atom topology from index_bonds
        # TODO: bachega's code
        initial = time.time()
        self._generate_topology_from_index_bonds()
        
        # Define molecules based on the atom topology
        # self.define_molecules()
        
        final = time.time()
        print('        Defining molecule indexes: ', final - initial)
        
        # Define Calpha backbone atoms
        self.define_Calpha_backbone()
    
    def find_bonded_and_nonbonded_atoms(self, selection, frame=0, gridsize=1.2,
                                        maxbond=2.4, tolerance=1.4):
        """
        Function doc: Determines bonded and nonbonded atoms based on selection in a VismolObject.
        
        Parameters:
        - selection: A dictionary containing the selected atoms (optional).
        - frame: Frame number for the atom coordinates (default is 0).
        - gridsize: Grid size for bond calculation (usually should not be changed, default is 0.8).
        - maxbond: Size of the maximum bond to be monitored (bonds greater than maxbond may be disregarded, default is 2.4 Å).
        - tolerance: Safety factor that multiplies the term (ra_cov + rb_cov)**2 (default is 1.4).
        - internal: If True, calculates the bindings for the object itself (used in the object's genesis).
                    If False, returns a list of atom_ids representing the bonds between atoms.
        
        Returns:
        If internal is True, the function updates the VismolObject's internal attributes:
        - index_bonds: List of pairs of atom indexes representing the bonded atoms.
        - bonds: List of Bond objects representing the bonds between atoms.
        - topology: Dictionary representing the atom topology and connectivity.
        - Non-bonded interactions and other attributes are also updated.
        
        If internal is False, the function returns a list of atom_ids representing the bonds between atoms.
        
        Note: This function calculates bonds between atoms based on their positions and covalent radii.
        
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
        # Initialize variables and data structures
        # initial = time.time()
        # atoms_frame_mask = np.zeros(len(self.atoms), bool)
        assert self.cov_radii_array is not None
            # self.cov_radii_array = np.empty(len(self.atoms), dtype=np.float32)
            # for i, atom in self.atoms.items():
            #     self.cov_radii_array[i] = atom.cov_rad
        
        # # Create a mask to identify atoms in the frame (all atoms if selection is None)
        # if selection is None:
        #     selection = self.atoms
        #     atoms_frame_mask[:] = True
        # else:
        #     atoms_frame_mask[:] = False
        #     for atom in selection.values():
        #         atoms_frame_mask[atom.atom_id] = True
        
        # Extract relevant data for bond calculation
        #cov_rads = self.cov_radii_array[atoms_frame_mask]
        #coords = self.frames[frame][atoms_frame_mask]
        cov_rads = self.cov_radii_array
        coords = self.frames[frame]
        
        # Generate indexes and grid positions for each atom in the selection
        indexes = []
        gridpos_list = []
        for atom in selection.values():
            indexes.append(atom.atom_id)
            gridpos_list.append(atom.get_grid_position(gridsize=gridsize, frame=frame))
        
        # logger.debug("Time used for preparing the atom mask, covalent radii list "\
        #                  "and grid positions: {}".format(time.time() - initial))
        
        # Calculate bonds based on grid positions and covalent radii and return the results
        index_bonds = cdist.get_atomic_bonds_from_grid(indexes, coords, cov_rads,
                                    gridpos_list, gridsize, maxbond, tolerance)
        return index_bonds
    
    def _bonds_from_pair_of_indexes_list(self, exclude_list = [['H','H']]):
        """ 
        Creates Bond objects based on pairs of indexes in self.index_bonds list.
        The bonds list is populated with the created Bond objects, and each
        atom involved in a bond is updated with the respective Bond object.
        
        self.index_bonds = [0,1  , 0,4  ,  1,3  , ...]
        self.bonds = [bond1(obj), bond2(obj), bond3(obj), ...] 
        """
        # assert self.bonds is None # Ensure the bonds list is not already initialized
        self.bonds = [] # Initialize an empty list to store the Bond objects
        
        new_index_bonds = []
        # Loop through the self.index_bonds list in pairs
        for i in range(0, len(self.index_bonds)-1, 2):
            
            index_i = self.index_bonds[i]    # Get the first atom's index of the bond
            index_j = self.index_bonds[i+1]  # Get the second atom's index of the bond
            
            is_excluded = False
            
            for excluded_bond in exclude_list: 
                if self.atoms[index_i].symbol in excluded_bond and self.atoms[index_j].symbol in excluded_bond:
                    is_excluded = True
            
            
            if is_excluded:
                pass
            else:
                new_index_bonds.append(index_i)
                new_index_bonds.append(index_j)
            
                #index_i = self.index_bonds[i]    # Get the first atom's index of the bond
                #index_j = self.index_bonds[i+1]  # Get the second atom's index of the bond
                
                
                
                # Create a Bond object with the atoms and their indexes
                bond = Bond(atom_i=self.atoms[index_i], atom_index_i=index_i,
                            atom_j=self.atoms[index_j], atom_index_j=index_j)
                
                # Add the created Bond object to the bonds list
                self.bonds.append(bond)
                
                # Update the atoms with the created Bond object, indicating their bond connections
                self.atoms[index_i].bonds.append(bond)
                self.atoms[index_j].bonds.append(bond)
        
         # Convert the index_bonds list to a numpy array of unsigned 32-bit integers
        self.index_bonds = new_index_bonds
        #print(self.index_bonds)
        self.index_bonds = np.array(self.index_bonds, dtype=np.uint32)
    
    def _get_non_bonded_from_bonded_list(self):
        """ Function doc """
        # assert self.non_bonded_atoms is None
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
                residue._is_protein()
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



