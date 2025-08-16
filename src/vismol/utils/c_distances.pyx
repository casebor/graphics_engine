#!/usr/bin/env python3
#cython: language_level=3
# -*- coding: utf-8 -*-


import multiprocessing
import numpy as np
cimport numpy as np


cpdef list calculate_grid_offset(gridsize, maxbond=2.6):
    '''
    This function calculates the offset for each grid element in the atomic grid. 
    The grid elements define regions in 3D space within a certain distance from 
    each atom. The gridsize parameter controls the size of each grid element, 
    and the maxbond parameter determines the maximum distance within which bonds 
    are considered. The function returns a list of grid offsets representing 
    grid elements.
    '''
    
    grid_offset_full = []
    borderGrid  = maxbond / gridsize
    borderGrid  = int(borderGrid)
    #------------------------- first floor -----------------------------
    #               |-------|-------|-------|-------|-------| 
    #               |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #               |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #               |-2,2,0 |-1,2,0 | 0,2,0 | 1,2,0 | 2,2,0 | 
    #               |-------|-------|-------|-------|-------| 
    #               |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #               |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #               |-1,1,0 |-1,1,0 | 0,1,0 | 1,1,0 | 2,1,0 | 
    #               |-------|-------|-------|-------|-------| 
    #               |       |       |XXXXXXX|\\\\\\\|\\\\\\\|
    #               |       |       |XXXXXXX|\\\\\\\|\\\\\\\|
    #               |-1,0,0 |-1,0,0 | 0,0,0 | 1,0,0 | 2,0,0 |
    #               |-------|-------|-------|-------|-------|
    #               |       |       |       |       |       |
    #               |       |       |       |       |       |
    #               |-1,-1,0|-1,-1,0| 0,-1,0| 1,-1,0| 2,-1,0|
    #               |-------|-------|-------|-------|-------|
    N = 0
    for i in range (-borderGrid, borderGrid + 1):
        n = 0
        for j in range(0, borderGrid + 1):
            #counter = i + n #+ 2
            
            # we don't need all the elements to the first floor
            if i < -1 and j == 0:
                pass
            
            else:
                # we don't need the [0,0,0] element to the first floor
                if [i,j, 0] == [0,0,0]:
                    pass
                else:
                    grid_offset_full.append([i,j,0])
                    N+=1
            
            n+=1
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
    
    # excluded_list = [] 
    # n = 0
    # for i in range (-borderGrid, borderGrid + 1):
    #     for j in range(-borderGrid, borderGrid + 1):
            
    #         for k in range(1, borderGrid + 1):
    #         #for k in range(0, borderGrid + 1):
    #             if [i, j,  k] in excluded_list:
    #                 pass
    #             else:
    #                 grid_offset_full.append([i,j,k])
    #                 #print([i,j,k], n+1)
    #             n+=1
    
    for i in range (-borderGrid, borderGrid + 1):
        for j in range(-borderGrid, borderGrid + 1):
            for k in range(1, borderGrid + 1):
                grid_offset_full.append([i,j,k])
    return grid_offset_full

cpdef double calculate_sqrt_distance(int i, int j, coords):
    """ 
    This Cython function calculates the squared distance between two atoms 
    with indices i and j using their coordinates coords. It is used to 
    calculate the distance between atoms and determine if a bond exists 
    between them. 
    """
    # dX              = (coords[i*3  ] - coords[j*3  ])**2
    # dY              = (coords[i*3+1] - coords[j*3+1])**2
    # dZ              = (coords[i*3+2] - coords[j*3+2])**2
    dX, dY, dZ = coords[i] - coords[j]
    dX = dX**2
    dY = dY**2
    dZ = dZ**2
    r_ij = dX + dY + dZ
    return r_ij

cpdef list get_connections_within_grid_element(list list_of_atoms, coords, cov_rad, double tolerance, gridsize):
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

    bonds_pair_of_indexes = []
    
    cdef double r_ij
    cdef int i
    cdef int atom_idx_i
    cdef int atom_idx_j
    cdef double cov_rad_ij_sqrt

    
    for i, atom_idx_i in enumerate(list_of_atoms[:-1]):
        for j, atom_idx_j in enumerate(list_of_atoms[i:]):    
            if atom_idx_i == atom_idx_j :
                pass            
            else:
                # print(atom_idx_i, "@", atom_idx_j, "@", coords)
                r_ij            = calculate_sqrt_distance(atom_idx_i, atom_idx_j, coords)
                # r_ij = np.linalg.norm(coords[atom_idx_i] - coords[atom_idx_j])
                cov_rad_ij_sqrt = ( (cov_rad[atom_idx_i] + cov_rad[atom_idx_j] )**2)*tolerance
                #cov_rad_ij_sqrt = ( (cov_rad[i] + cov_rad[j] )**2)*1.4
                
                #print (atom_idx_i, atom_idx_j,cov_rad[atom_idx_j],cov_rad[atom_idx_j] , r_ij, cov_rad_ij_sqrt)
                if r_ij <= cov_rad_ij_sqrt:
                    pass
                    #bonds_pair_of_indexes.append([atom_idx_i , atom_idx_j])
                    bonds_pair_of_indexes.append(atom_idx_i)
                    bonds_pair_of_indexes.append(atom_idx_j)
                else:
                    pass

    return bonds_pair_of_indexes


cpdef list get_connections_between_grid_elements(list atomic_grid1, 
                                                       list atomic_grid2, 
                                                       coords, 
                                                       cov_rad, 
                                                       double  tolerance, 
                                                       gridsize):
    '''
    This function calculates the bonds between atoms in two different grid
    elements. It takes two atomic grids, their coordinates, covalent radii, 
    a tolerance factor, and the grid size. It returns a list of pairs of atom
    indices representing the bonds between the two grid elements.
    '''
    cdef double r_ij
    cdef double cov_rad_ij_sqrt
    cdef int atom_idx_i
    cdef int atom_idx_j

    cdef list bonds_pair_of_indexes
    bonds_pair_of_indexes = []
    
    
    if atomic_grid1 == atomic_grid2:
        pass
    else:
        for atom_idx_i in atomic_grid1:   
            #xyz       = atom_i[2]
            #atom_ix   = xyz[0]
            #atom_iy   = xyz[1]    
            #atom_iz   = xyz[2]
            #cov_rad_i = atom_i[3]
            #index_i   = atom_i[0]

            for atom_idx_j in atomic_grid2:    

                if atom_idx_i == atom_idx_j :
                    pass            
                else:
                    r_ij            = calculate_sqrt_distance(atom_idx_i, atom_idx_j, coords)
                    # r_ij = np.linalg.norm(coords[atom_idx_i] - coords[atom_idx_j])
                    #cov_rad_ij_sqrt = ( (cov_rad[atom_idx_i] + cov_rad[atom_idx_j] )**2)*1.4
                    cov_rad_ij_sqrt = ( (cov_rad[atom_idx_i] + cov_rad[atom_idx_j] )**2)*tolerance
                    #print (atom_idx_i, atom_idx_j,cov_rad[atom_idx_j],cov_rad[atom_idx_j] , r_ij, cov_rad_ij_sqrt)
                    if r_ij <= cov_rad_ij_sqrt:
                        #bonds_pair_of_indexes.append([atom_idx_i , atom_idx_j])
                        bonds_pair_of_indexes.append(atom_idx_i )
                        bonds_pair_of_indexes.append(atom_idx_j )
                    else:
                        pass

    return bonds_pair_of_indexes


cpdef dict build_the_atomic_grid ( list indexes     ,
                                         list gridpos_list):
    
    '''
    This function builds the atomic grid by grouping atoms based on their 
    grid positions. It takes a list of atom indices and their corresponding 
    grid positions. It returns a dictionary where each key represents a grid 
    element, and the value is a list of atom indices within that grid element.
    '''
    
    cdef int atom

    atomic_grid = {}
    
    for atom, grid_pos in enumerate(gridpos_list):
        if grid_pos in atomic_grid:
            atomic_grid[grid_pos].append(indexes[atom])
        else:
            
            atomic_grid[grid_pos] = []
            atomic_grid[grid_pos].append(indexes[atom])
    
    return atomic_grid


cpdef list get_atomic_bonds_from_grid(list indexes, coords, cov_rad, list gridpos_list,
                                      double gridsize, double maxbond, double tolerance):
    '''
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
    '''
    
    #cdef double tolerance
    ##cpdef double maxbond
    ##maxbond   = 3.0
    #tolerance = 1.8
    
    atomic_grid = build_the_atomic_grid ( indexes     ,
                                                gridpos_list)

    grid_offset_full = calculate_grid_offset(gridsize, maxbond)
   
    #print ('gridsize',  gridsize )
    #print ('maxbond',   maxbond )
    #print ('borderGrid',  int(maxbond/gridsize)   )
    #print ('grid_offset_full', len(grid_offset_full))
    
    n = 0
    bonds_pair_of_indexes = [] 
   
    for element in atomic_grid.keys():
        
        atomic_grid1   = atomic_grid[element]

        #----------------------------------------------------------------------------------------#
        if len(atomic_grid1) == 1:
            pass
        else:
            pass
            element1_bonds = get_connections_within_grid_element( atomic_grid1, coords, cov_rad, tolerance, gridsize)        
            bonds_pair_of_indexes += element1_bonds
        #----------------------------------------------------------------------------------------#
        n+=1
        #----------------------------------------------------------------------------------------#
        for offset_element in  grid_offset_full:              
            
            element2  = (element[0]+offset_element[0], 
                         element[1]+offset_element[1], 
                         element[2]+offset_element[2]) 

            if element2 in atomic_grid:                        
                n+=1
                pass
                element1_2_bonds = get_connections_between_grid_elements(atomic_grid1, atomic_grid[element2], coords, cov_rad, tolerance, gridsize)
                bonds_pair_of_indexes += element1_2_bonds
        #----------------------------------------------------------------------------------------#
    #print('n interaction', n)
    return bonds_pair_of_indexes
