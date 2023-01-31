#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import multiprocessing
import numpy as np
cimport numpy as np


cpdef list calculate_grid_offset(gridsize, maxbond=2.6):
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
    
    excluded_list = [] 
    n = 0
    for i in range (-borderGrid, borderGrid + 1):
        for j in range(-borderGrid, borderGrid + 1):
            
            for k in range(1, borderGrid + 1):
            #for k in range(0, borderGrid + 1):
                if [i, j,  k] in excluded_list:
                    pass
                else:
                    grid_offset_full.append([i,j,k])
                    #print([i,j,k], n+1)
                n+=1
    return grid_offset_full

cpdef double calculate_sqrt_distance(int i, int j, coords):
    """ Function doc """
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
    
    cpdef double r_ij
    cpdef int i
    cpdef int atom_idx_i
    cpdef int atom_idx_j
    cpdef double cov_rad_ij_sqrt

    
    for i, atom_idx_i in enumerate(list_of_atoms[:-1]):
        for j, atom_idx_j in enumerate(list_of_atoms[i:]):    
            if atom_idx_i == atom_idx_j :
                pass            
            else:
                # print(atom_idx_i, "@", atom_idx_j, "@", coords)
                r_ij            = calculate_sqrt_distance(atom_idx_i, atom_idx_j, coords)
                # r_ij = np.linalg.norm(coords[atom_idx_i] - coords[atom_idx_j])
                cov_rad_ij_sqrt = ( (cov_rad[atom_idx_i] + cov_rad[atom_idx_j] )**2)*1.4
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
    
    cpdef double r_ij
    cpdef double cov_rad_ij_sqrt
    cpdef int atom_idx_i
    cpdef int atom_idx_j

    cpdef list bonds_pair_of_indexes
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
                    cov_rad_ij_sqrt = ( (cov_rad[atom_idx_i] + cov_rad[atom_idx_j] )**2)*1.4
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
    cpdef int atom

    atomic_grid = {}
    
    for atom, grid_pos in enumerate(gridpos_list):
        if grid_pos in atomic_grid:
            atomic_grid[grid_pos].append(indexes[atom])
        else:
            
            atomic_grid[grid_pos] = []
            atomic_grid[grid_pos].append(indexes[atom])
    
    return atomic_grid


cpdef list get_atomic_bonds_from_grid(list indexes, coords, cov_rad, list gridpos_list,
                                      double gridsize, double maxbond):
    
    
    cpdef double tolerance
    #cpdef double maxbond
    #maxbond   = 3.0
    tolerance = 1.4
    
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
