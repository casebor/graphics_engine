#!/usr/bin/env python3
#cython: language_level=3
# -*- coding: utf-8 -*-


import numpy as np
cimport numpy as np
from vismol.model.bond import Bond


cpdef list get_bond_pairs(molecule, np.ndarray cov_radii_array,
                     np.ndarray atom_ids, np.float32_t gridsize=1.2,
                     np.float32_t maxbond=2.4, np.float32_t tolerance=1.4):
    """
    function is responsible for calculating the bonds between atoms within a given grid. It takes the following parameters:

    indexes: A list of atom indices.
    coords:  A list of atom coordinates.
    cov_rad: A list of covalent radii for each atom.
    
    gridpos_list: A list of grid positions that represent the atomic grid.
    gridsize:     A double representing the grid size.
    maxbond:      A double representing the maximum bond length.
    
    
    The function first builds an atomic grid using the build_the_atomic_grid 
    function, which organizes atoms into grid elements based on their positions. 
    It then calculates the grid offset using the calculate_grid_offset function. 
    The grid offset is used to define neighboring grid elements around each 
    element in the atomic grid.

    The function iterates over each grid element and calculates the bonds 
    between atoms within that element. If there is more than one atom in 
    the element, it calls the get_connections_within_grid_element function 
    to calculate the intra-element bonds. The result is added to the 
    bonds_pair_of_indexes list.

    Next, the function iterates over neighboring grid elements based on the 
    grid offset and calculates the bonds between atoms in the current element 
    and those in the neighboring element. The result is again added to the 
    bonds_pair_of_indexes list.

    Finally, the function returns the bonds_pair_of_indexes list containing 
    the pairs of atom indices representing the bonds within the atomic grid.
    """
    assert gridsize > 0.0, "Value error: gridsize must be larger than 0.0"
    cdef list gridpos_list = []
    cdef list new_atom_ids = []
    cdef np.ndarray coords = molecule.frames[0]
    
    for i in range(coords.shape[0]):
        if molecule.atoms[i].altloc is not None:
            if molecule.atoms[i].altloc == "A":
                gridpos_list.append(get_grid_position(coords[i], gridsize))
                new_atom_ids.append(i)
        else:
            gridpos_list.append(get_grid_position(coords[i], gridsize))
            new_atom_ids.append(i)
    
    cdef dict atomic_grid = build_atomic_grid(np.array(new_atom_ids, dtype=np.uint32), gridpos_list)
    cdef list grid_offset_full = calculate_grid_offset(gridsize, maxbond)
    
    cdef list atomic_grid1, element1_bonds, element1_2_bonds
    cdef tuple atom_grid, offset_grid
    cdef list bonds_pair_of_indexes = []
    for atom_grid in atomic_grid.keys():
        atomic_grid1 = atomic_grid[atom_grid]
        if len(atomic_grid1) == 1:
            pass
        else:
            element1_bonds = get_connections_within_grid_element(atomic_grid1,
                                        coords, cov_radii_array, tolerance)
            bonds_pair_of_indexes.extend(element1_bonds)
        for offset_grid in grid_offset_full:
            element2  = (atom_grid[0]+offset_grid[0],
                         atom_grid[1]+offset_grid[1],
                         atom_grid[2]+offset_grid[2])
            if element2 in atomic_grid:
                element1_2_bonds = get_connections_between_grid_elements(atomic_grid1,
                        atomic_grid[element2], coords, cov_radii_array, tolerance)
                bonds_pair_of_indexes.extend(element1_2_bonds)
    return bonds_pair_of_indexes


cpdef tuple get_grid_position(np.ndarray coords, np.float32_t gridsize):
    """ Function doc """
    assert gridsize > 0.0, "Value error: gridsize must be larger than 0.0"
    cdef tuple gridpos = (int(coords[0]/gridsize), int(coords[1]/gridsize), int(coords[2]/gridsize))
    return gridpos


cpdef dict build_atomic_grid(np.ndarray atom_ids, list gridpos_list):
    """
    This function builds the atomic grid by grouping atoms based on their 
    grid positions. It takes a list of atom indices and their corresponding 
    grid positions. It returns a dictionary where each key represents a grid 
    element, and the value is a list of atom indices within that grid element.
    """
    cdef np.int32_t atom_i
    cdef dict atomic_grid = {}
    for atom_i, grid_pos in enumerate(gridpos_list):
        if grid_pos in atomic_grid:
            atomic_grid[grid_pos].append(atom_ids[atom_i])
        else:
            atomic_grid[grid_pos] = []
            atomic_grid[grid_pos].append(atom_ids[atom_i])
    return atomic_grid


cpdef list calculate_grid_offset(np.float32_t gridsize, np.float32_t maxbond):
    """
    This function calculates the offset for each grid element in the atomic grid. 
    The grid elements define regions in 3D space within a certain distance from 
    each atom. The gridsize parameter controls the size of each grid element, 
    and the maxbond parameter determines the maximum distance within which bonds 
    are considered. The function returns a list of grid offsets representing 
    grid elements.
    """
    assert gridsize > 0.0, "Value error: gridsize must be larger than 0.0"
    cdef list grid_offset_full = []
    cdef np.int32_t border_grid, i, j, k
    border_grid = int(maxbond/gridsize)
    #---------------- first floor -------------
    # |-------|-------|-------|-------|-------|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |-2,2,0 |-1,2,0 | 0,2,0 | 1,2,0 | 2,2,0 |
    # |-------|-------|-------|-------|-------|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |-1,1,0 |-1,1,0 | 0,1,0 | 1,1,0 | 2,1,0 |
    # |-------|-------|-------|-------|-------|
    # |       |       |XXXXXXX|\\\\\\\|\\\\\\\|
    # |       |       |XXXXXXX|\\\\\\\|\\\\\\\|
    # |-1,0,0 |-1,0,0 | 0,0,0 | 1,0,0 | 2,0,0 |
    # |-------|-------|-------|-------|-------|
    # |       |       |       |       |       |
    # |       |       |       |       |       |
    # |-1,-1,0|-1,-1,0| 0,-1,0| 1,-1,0| 2,-1,0|
    # |-------|-------|-------|-------|-------|
    
    for i in range(-border_grid, border_grid+1):
        for j in range(0, border_grid+1):
            # we don't need all the elements to the first floor
            if i < -1 and j == 0:
                pass
            else:
                # we don't need the [0,0,0] element to the first floor
                if [i, j, 0] == [0, 0, 0]:
                    pass
                else:
                    grid_offset_full.append((i, j, 0))
    
    #--------------------- floors above---------------------------------
    #                                     |-------|-------|-------|-------|
    #                                     |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    #                                     |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    #                                     |-1,2,2 | 0,2,2 | 1,2,2 | 2,2,2 |
    #                                     |-------|-------|-------|-------|
    # |-------|-------|-------|-------|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |-1,1,2 | 0,1,2 | 1,1,2 | 2,1,2 |
    # |-1,2,1 | 0,2,1 | 1,2,1 | 2,2,1 |   |-------|-------|-------|-------|
    # |-------|-------|-------|-------|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |-1,0,2 | 0,0,2 | 1,0,2 | 2,0,2 |
    # |-1,1,1 | 0,1,1 | 1,1,1 | 2,1,1 |   |-------|-------|-------|-------|
    # |-------|-------|-------|-------|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |-1,-1,2| 0,-1,2| 1,-1,2| 2,-1,2|
    # |-1,0,1 | 0,0,1 | 1,0,1 | 2,0,1 |   |-------|-------|-------|-------|
    # |-------|-------|-------|-------|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |-1,-1,1| 0,-1,1| 1,-1,1| 2,-1,1|
    # |-------|-------|-------|-------|
    
    # TODO: Is this part used anywhere? The variable excluded_list is always empty
    #       and the tuples are always added to the grid_offset_full list
    for i in range (-border_grid, border_grid+1):
        for j in range(-border_grid, border_grid+1):
            for k in range(1, border_grid+1):
                grid_offset_full.append((i, j, k))
    return grid_offset_full


cpdef list get_connections_within_grid_element(list atom_grid, np.ndarray coords,
                np.ndarray cov_rads, np.float32_t tolerance):
    """
    This function calculates the bonds between atoms within a single 
    element of the atomic grid. It takes a list of atoms, their coordinates, 
    covalent radii (cov_rad), a tolerance factor, and the grid size. 
    It returns a list of pairs of atom indices representing the bonds 
    within the grid element.
    
        Calculate the distances and bonds 
        between atoms within a single element 
        of the atomic grid
        
                  |-------|-------|-------|
                  |       |       |       |
                  |       |       |       |
                  |       |       |       |
                  |-------|-atoms-|-------|
                  |       |       |       |
                  |       | i<->j |       |
                  |       |       |       |
                  |-------|-------|-------|
                  |       |       |       |
                  |       |       |       |
                  |       |       |       |
                  |-------|-------|-------|
    
    
    
    atoms = [[index, at_name, cov_rad,  at_pos, at_res_i, at_res_n, at_ch], ...]
            each elemte is a list contain required data.
    
    
    bonds_pair_of_indexes [[a,b],[b,c], ...] where a and b are indices. 
    returns a list of pair of indices "bonds_pair_of_indexes"
    """
    cdef list bonds_pair_of_indexes = []
    cdef np.int32_t i, j, atom_idx_i, atom_idx_j
    cdef np.float32_t r_ij, cov_rad_ij_sqrd
    for i, atom_idx_i in enumerate(atom_grid[:-1]):
        for j, atom_idx_j in enumerate(atom_grid[i:]):
            if atom_idx_i != atom_idx_j:
                r_ij = get_sqrd_distance(atom_idx_i, atom_idx_j, coords)
                cov_rad_ij_sqrd = ((cov_rads[atom_idx_i] + cov_rads[atom_idx_j])**2)*tolerance
                if r_ij <= cov_rad_ij_sqrd:
                    bonds_pair_of_indexes.append(atom_idx_i)
                    bonds_pair_of_indexes.append(atom_idx_j)
    return bonds_pair_of_indexes


cpdef np.float32_t get_sqrd_distance(np.int32_t i, np.int32_t j, np.ndarray coords):
    """
    This Cython function calculates the squared distance between two atoms 
    with indices i and j using their coordinates coords. It is used to 
    calculate the distance between atoms and determine if a bond exists 
    between them.
    """
    cdef np.float32_t dx, dy, dz, r_ij
    dx, dy, dz = coords[i] - coords[j]
    dx = dx**2
    dy = dy**2
    dz = dz**2
    r_ij = dx + dy + dz
    return r_ij


cpdef list get_connections_between_grid_elements(list atomic_grid1, list atomic_grid2,
                np.ndarray coords, np.ndarray cov_rads, np.float32_t tolerance):
    """
    This function calculates the bonds between atoms in two different grid
    elements. It takes two atomic grids, their coordinates, covalent radii, 
    a tolerance factor, and the grid size. It returns a list of pairs of atom
    indices representing the bonds between the two grid elements.
    """
    cdef np.float32_t r_ij, cov_rad_ij_sqrd
    cdef np.int32_t atom_idx_i, atom_idx_j
    cdef list bonds_pair_of_indexes = []
    
    if atomic_grid1 != atomic_grid2:
        for atom_idx_i in atomic_grid1:
            for atom_idx_j in atomic_grid2:
                if atom_idx_i != atom_idx_j:
                    r_ij = get_sqrd_distance(atom_idx_i, atom_idx_j, coords)
                    cov_rad_ij_sqrd = ((cov_rads[atom_idx_i] + cov_rads[atom_idx_j])**2)*tolerance
                    if r_ij <= cov_rad_ij_sqrd:
                        bonds_pair_of_indexes.append(atom_idx_i)
                        bonds_pair_of_indexes.append(atom_idx_j)
    return bonds_pair_of_indexes


cpdef np.ndarray bonds_from_pair_of_indexes_list(molecule, list index_bonds, list exclude_list=None):
    """
    Creates Bond objects based on pairs of indexes in self.index_bonds list.
    The bonds list is populated with the created Bond objects, and each
    atom involved in a bond is updated with the respective Bond object.
    
    self.index_bonds = [0,1  , 0,4  ,  1,3  , ...]
    self.bonds = [bond1(obj), bond2(obj), bond3(obj), ...] 
    """
    if exclude_list is None:
        exclude_list = [["H","H"]]
    
    cdef set bonds = set() # Initialize an empty list to store the Bond objects
    cdef list new_index_bonds = []
    cdef np.int32_t i, index_i, index_j
    cdef bint is_excluded
    # Loop through the self.index_bonds list in pairs
    for i in range(0, len(index_bonds)-1, 2):
        
        index_i = index_bonds[i]    # Get the first atom's index of the bond
        index_j = index_bonds[i+1]  # Get the second atom's index of the bond
        is_excluded = False
        
        for excluded_bond in exclude_list:
            if molecule.atoms[index_i].symbol in excluded_bond and \
               molecule.atoms[index_j].symbol in excluded_bond:
                is_excluded = True
        
        if not is_excluded:
            # Create a Bond object with the atoms and their indexes
            bond = Bond(atom_i=molecule.atoms[index_i],
                        atom_j=molecule.atoms[index_j])
            # Add the created Bond object to the bonds list
            bonds.add(bond)
            # Update the atoms with the created Bond object, indicating their bond connections
            molecule.atoms[index_i].bonds.add(bond)
            molecule.atoms[index_j].bonds.add(bond)
    if molecule.has_altloc:
        alt_bonds = copy_altloc_topology(molecule)
        bonds.update(alt_bonds)
    new_index_bonds.extend(_get_indexes_from_bonds(bonds))
    molecule.topology.bonds = bonds
     # Convert the index_bonds list to a numpy array of unsigned 32-bit integers
    return np.array(new_index_bonds, dtype=np.uint32)


cpdef np.ndarray get_non_bonded_from_bonded_list(molecule, np.ndarray index_bonds):
    """ Function doc """
    bonded_set = set(index_bonds)
    cdef list non_bonded_atoms = []
    for i, atom in molecule.atoms.items():
        if i in bonded_set:
            atom.nonbonded = False
        else:
            atom.nonbonded = True
            non_bonded_atoms.append(i)
    non_bonded_atoms.sort()
    return np.array(non_bonded_atoms, dtype=np.uint32)


cpdef list copy_altloc_topology(molecule):
    cdef list alt_bonds = []
    cdef dict altlocs
    for chain in molecule.chains.values():
        for res in chain.residues.values():
            if res.has_altloc:
                altlocs = {}
                for atom in res.atoms.values():
                    if atom.altloc in altlocs:
                        altlocs[atom.altloc].append(atom)
                    else:
                        altlocs[atom.altloc] = [atom]
                pattern = altlocs["A"]
                del altlocs["A"]
                for alt_atoms in altlocs.values():
                    for pat_at in pattern:
                        new_bonds, new_indexes = _copy_atom_bonds(pat_at, alt_atoms)
                        alt_bonds.extend(new_bonds)
    return alt_bonds


cpdef list _get_indexes_from_bonds(set bonds):
    cdef list indexes = []
    for bond in bonds:
        indexes.append(bond.atom_index_i)
        indexes.append(bond.atom_index_j)
    return indexes


cpdef tuple _copy_atom_bonds(pat_at, list alt_atoms):
    cdef set bonds = set()
    cdef list indexes = []
    for bond in pat_at.bonds:
        if bond.atom_i.name == pat_at.name:
            pos_i = _get_name_pos(bond.atom_i, alt_atoms)
            pos_j = _get_name_pos(bond.atom_j, alt_atoms)
            alt_i = alt_atoms[pos_i]
            if pos_j == -1:
                alt_j = bond.atom_j
            else:
                alt_j = alt_atoms[pos_j]
        elif bond.atom_j.name == pat_at.name:
            pos_i = _get_name_pos(bond.atom_i, alt_atoms)
            pos_j = _get_name_pos(bond.atom_j, alt_atoms)
            alt_j = alt_atoms[pos_j]
            if pos_i == -1:
                alt_i = bond.atom_i
            else:
                alt_i = alt_atoms[pos_i]
        else:
            pass
        new_bond = Bond(atom_i=alt_i, atom_j=alt_j)
        bonds.add(new_bond)
        indexes.append(alt_i.atom_id)
        indexes.append(alt_j.atom_id)
        alt_i.bonds.add(new_bond)
        alt_j.bonds.add(new_bond)
    return (bonds, indexes)


cdef np.int32_t _get_name_pos(pat_atom, list atoms):
    for i, atom in enumerate(atoms):
        if pat_atom.residue != atom.residue:
            return -1
        if atom.name == pat_atom.name:
            return i
    return -1

