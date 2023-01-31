#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import multiprocessing
import numpy as np
cimport numpy as np



cpdef tuple selection_spherical_expansion (set selected_indexes, set selectable_indexes,  coordinates, float radius ):
        cpdef float radius_sqr
        cpdef float dx
        cpdef float dy
        cpdef float dz
        
        radius_sqr = radius**2
        new_selected_indexes   = set()
        
        for i in selected_indexes:
            #print ('len(selectable_indexes)', len(selectable_indexes))
            for j in selectable_indexes:
                
                i_xyz  = coordinates[i]
                j_xyz  = coordinates[j]
                
                dx = (i_xyz[0] - j_xyz[0])**2
                dy = (i_xyz[1] - j_xyz[1])**2
                dz = (i_xyz[2] - j_xyz[2])**2
                
                if dx + dy + dz <= radius_sqr:
                    new_selected_indexes.add(j)
                else:
                    pass
            
            selectable_indexes = selectable_indexes - new_selected_indexes
        
        for index in selected_indexes:
            new_selected_indexes.add(index)

        return new_selected_indexes,  selectable_indexes



