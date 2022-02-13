#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  molecular_model.py
#  
#  Copyright 2016 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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

import numpy as np
import vModel.Vectors as LA
from vModel import VismolObject

class VisMolPickingSelection:
    """ Class doc """
    
    def __init__ (self, vismolSession):
        """ Class initialiser """
        self.picking_selections_list = [None]*4
        self.picking_selections_list_index = []
        self.vismolSession = vismolSession

    
    def _generate_picking_selection_coordinates (self):
        """ Function doc """
        pass
        #for i,atom in enumerate(self.picking_selections_list):
        #    if atom is not None:
        #        coord = [atom.Vobject.frames[frame][(atom.index-1)*3  ],
        #                 atom.Vobject.frames[frame][(atom.index-1)*3+1],
        #                 atom.Vobject.frames[frame][(atom.index-1)*3+2],]
        #                
        #        rep.draw_selected(atom, coord, [0.83, 0.48, 1])
        #        rep.draw_numbers(atom, i+1, coord)
    
    
    
    
    def selection_function_picking (self, selected):
        """ Function doc """
        if selected is None:
            self.picking_selections_list = [None]*len(self.picking_selections_list)
            #self.selected_residues = []
        else:
            if selected not in self.picking_selections_list:
                for i in range(len(self.picking_selections_list)):
                    if self.picking_selections_list[i] == None:
                        self.picking_selections_list[i] = selected
                        selected = None
                        break
                if selected is not None:
                    self.picking_selections_list[len(self.picking_selections_list)-1] = selected
            else:
                for i in range(len(self.picking_selections_list)):
                    if self.picking_selections_list[i] == selected:
                        self.picking_selections_list[i] = None

        
        #print('\n\nDistances:')
        c = 0
        for atom1 in self.picking_selections_list:
            for atom2 in self.picking_selections_list[c+1:]:
            
                if atom1 and atom2:
                    dist = self.vismolSession._get_distance_atom1_atom2 ( atom1, atom2 )
                    name1 = atom1.name
                    name2 = atom2.name
                    print ('atom',name1, 'atom',name2,  dist)
                    
            
            c += 1
        
        atom1 = self.picking_selections_list[0]
        atom2 = self.picking_selections_list[1]
        atom3 = self.picking_selections_list[2]
        atom4 = self.picking_selections_list[3]
        #print(atom1,atom2,atom3,atom4)
        #self.refresh_pk1pk2_representations( vobj_label ='pk1pk2', atom1 = atom1, atom2 = atom2)
        
        if atom1 and atom2:
            #print ('line 95')
            self.refresh_pk1pk2_representations( vobj_label =  'pk1pk2',
                                                      atom1 = atom1    , 
                                                      atom2 = atom2    )
            self.vismolSession.vismol_geometric_object_dic['pk1pk2'].representations['dotted_lines'].active =  True
            if atom3:
                xyz1 = atom1.coords()
                xyz2 = atom2.coords()
                xyz3 = atom3.coords()
                
                xyz1 = [ xyz1[0] - xyz2[0], xyz1[1] - xyz2[1],   xyz1[2] - xyz2[2]]
                xyz3 = [ xyz3[0] - xyz2[0], xyz3[1] - xyz2[1],   xyz3[2] - xyz2[2]]

                angle = LA.angle(xyz1, xyz3)
                print ('Angle: ', angle*57.297)
                text =  'Angle: '+ str( angle*57.297)
                self.vismolSession.main_session.statusbar_main.push(1,text)
                if atom4:
                    xyz4 = atom4.coords()
                    angle = LA.dihedral(xyz1, xyz2, xyz3, xyz4)
                    print ('Dihedral: ', angle*57.297)
        
        else:
            ##print(self.vismolSession.vismol_geometric_object_dic['pk1pk2'],self.vismolSession.vismol_geometric_object_dic['pk1pk2'].active )
            if self.vismolSession.vismol_geometric_object_dic['pk1pk2']:
                #print('120')
                self.vismolSession.vismol_geometric_object_dic['pk1pk2'].representations['dotted_lines'].active =  False
                                                                            
                                                                            
        if atom2 and atom3:                                                 
            #print ('line 95')                                              
            self.refresh_pk1pk2_representations( vobj_label = 'pk2pk3' ,    
                                                      atom1 = atom2    ,    
                                                      atom2 = atom3    )    
            self.vismolSession.vismol_geometric_object_dic['pk2pk3'].representations['dotted_lines'].active =  True
        else:
            #print('128')
            if self.vismolSession.vismol_geometric_object_dic['pk2pk3']:
                self.vismolSession.vismol_geometric_object_dic['pk2pk3'].representations['dotted_lines'].active =  False
                                                                            
                                                                            
        if atom3 and atom4:                                                 
            self.refresh_pk1pk2_representations( vobj_label =  'pk3pk4',    
                                                      atom1 = atom3    ,    
                                                      atom2 = atom4    )    
            self.vismolSession.vismol_geometric_object_dic['pk3pk4'].representations['dotted_lines'].active =  True

        else:
            if self.vismolSession.vismol_geometric_object_dic['pk3pk4']:
                self.vismolSession.vismol_geometric_object_dic['pk3pk4'].representations['dotted_lines'].active =  False

    
    def refresh_pk1pk2_representations( self, vobj_label ='pk1pk2',
                                               atom1 = None   , 
                                               atom2 = None   ):
        #print('bulding vobject_picking line 121')

        xyz1 = atom1.coords()
        xyz2 = atom2.coords()
        frame = np.array(xyz1 + xyz2, dtype=np.float32)
        
        if self.vismolSession.vismol_geometric_object_dic[vobj_label]:
            #print('bulding vobject_picking line 125')
            self.vismolSession.vismol_geometric_object_dic[vobj_label].frames = [frame]
            #print('bulding vobject_picking line 133')
            self.vismolSession.vismol_geometric_object_dic[vobj_label].representations['dotted_lines']._make_gl_vao_and_vbos()
            self.vismolSession.vismol_geometric_object_dic[vobj_label].active = True
        else:
            #print('bulding vobject_picking line 131')
            atoms = []
            atoms.append({
                          'index'      : 0             , 
                          'name'       : 'pK'           , 
                          'resi'       : ''            , 
                          'resn'       : ''            , 
                          'chain'      : ''            , 
                          'symbol'     : 'pK'           , 
                          'occupancy'  : 00.00         , 
                          'bfactor'    : 0.00          , 
                          'charge'     : 0.00           
                          })

            atoms.append({
                          'index'      : 1     , 
                          'name'       : 'pK'   , 
                          'resi'       : ''    , 
                          'resn'       : ''    , 
                          'chain'      : ''    , 
                          'symbol'     : 'pK'   , 
                          'occupancy'  : 00.00 , 
                          'bfactor'    : 0.00  , 
                          'charge'     : 0.00   
                          })
            
            frame = np.array(xyz1 + xyz2, dtype=np.float32)
            #print('bulding vobject_picking line 161')
            self.vobject_picking = VismolObject.VismolObject(name                           = 'UNK'              , 
                                                             atoms                          = atoms              ,
                                                             vismolSession                  = self.vismolSession , 
                                                             trajectory                     = [frame],
                                                             bonds_pair_of_indexes          = [0,1] , 
                                                             auto_find_bonded_and_nonbonded = False)
            
            self.vobject_picking.active = True
            self.vobject_picking.set_model_matrix(self.vismolSession.glwidget.vm_widget.model_mat)
            
            self.vobject_picking.create_new_representation(rtype = 'dotted_lines')
            self.vismolSession.vismol_geometric_object_dic[vobj_label] = self.vobject_picking
            #print('bulding vobject_picking line 174')

        
        


class VisMolViewingSelection:
    """ Class doc """
    
    def __init__ (self, vismolSession):
        #---------------------------------------------------------------
        #                S E L E C T I O N S
        #---------------------------------------------------------------
        self.active = False
        
        self._selection_mode       = 'residue'
        self.selected_objects      = {} #dic of VisMol objects (obj)
        self.selected_atoms        = [] #List of atoms objects (obj)
        self.selected_atoms_coords = [] #coordinate (floats) x y z
        self.vismolSession             = vismolSession
    
    def get_selection_info (self):
        """ Function doc """
        ##print('self._selection_mode          ',self._selection_mode           )
        ##print('self.selected_objects         ',self.selected_objects          )
        #
        ##print('self.selected_atoms       ',self.selected_atoms        )
        #print('Selection defiend with  ',len(self.selected_atoms), 'atom(s)')
        ##print('self.selected_atoms_coords',self.selected_atoms_coords )
        
        '''
        for atom in self.selected_atoms:
            #print(atom.name,                             # nome
                  atom.resi,                             # residue index
                  atom.Vobject.residues[atom.resi].resn, # residue name  taken from residues dic
                  atom.residue)                          # residue obj
        '''
            
        ##print(atom.name, atom.resi, atom.Vobject.residues) 
        ##print('self.selected_residues        ',self.selected_residues         )
        ##print('self.selected_atoms           ',self.selected_atoms            )
        ##print('self.selected_frames          ',self.selected_frames           )
        return len(self.selected_atoms)
    
    def selecting_by_atom (self, selected_atom, disable = True):
        """
        The "disable" variable does not allow, if the selected 
        atom is already in the selected list, to be removed. 

        The disable variable is "False" for when we use 
        selection by area (selection box)
        """
        
        self._clear_selection_buffer(selected_atom)

        if selected_atom not in self.selected_atoms:
            self.selected_atoms.append(selected_atom)
            selected_atom.selected = True
            
        else:
            if disable:
                index = self.selected_atoms.index(selected_atom)
                self.selected_atoms.pop(index)
                selected_atom.selected = False
            else:
                pass
                
    def selecting_by_residue (self, selected_atom, disable = True):
        """       
        when disable = True
        If the object selection is disabled, all atoms within the 
        residue will be set to False 
        
        The disable variable does not allow, if the selected 
        atom is already in the selected list, to be removed. 
        
        The disable variable is "False" for when we use 
        selection by area (selection box)  """
        
        self._clear_selection_buffer(selected_atom)

            
        #if selected
        
        # if the selected atoms is not in the selected list
        if selected_atom not in self.selected_atoms:
            
            for atom in selected_atom.residue.atoms:
                '''print (len(selected.residue.atoms), atom.name, atom.index)'''
                
                # the atom is not on the list -  add atom by atom
                if atom not in self.selected_atoms:
                    self.selected_atoms.append(atom)
                    atom.selected = True
                # the atom IS on the list - do nothing 
                else:
                    pass
    
        # else: if the selected atoms IS in the selected list
        else:
            if disable:
                # So, add all atoms  - selected residue <- selected.resi
                for atom in selected_atom.residue.atoms:
                    
                    # the atom is not on the list -  add atom by atom
                    if atom in self.selected_atoms:
                        index = self.selected_atoms.index(atom)
                        self.selected_atoms.pop(index)
                        atom.selected = False

                    # the atom IS in the list - do nothing 
                    else:
                        pass   

            else:
                pass
                
    def selecting_by_chain (self, selected_atom, disable = True):
        ''' 
        when disable = True
        If the object selection is disabled, all atoms in the system will be set to False 
        
        The disable variable does not allow, if the selected 
        atom is already in the selected list, to be removed. 
        
        The disable variable is "False" for when we use 
        selection by area (selection box) 
        
        '''
        #vismolSession      = selected_atom.Vobject.vismolSession
        self._clear_selection_buffer(selected_atom)
        
        '''
        vismolSession      = selected_atom.Vobject.vismolSession
        
        #------------------------------------------------
        # Clearing the selections buffer
        # assigning false to all atoms loaded into memory
        if self.active:
            pass
        else:
            for Vobject in vismolSession.vismol_objects:
                for atom in Vobject.atoms:
                    atom.selected = False
        #------------------------------------------------
        '''
                    
        print ('selecting_by_chain', selected_atom )
        # if the selected atoms is not on the selected list
        
        if selected_atom.selected == False:
        
            chain = selected_atom.chain
            print (chain)
            
            for atom in selected_atom.Vobject.atoms_by_chains[chain]:
                atom.selected = True 
            
        else:
            if disable:
                chain = selected_atom.chain
                for atom in selected_atom.Vobject.atoms_by_chains[chain]:
                    atom.selected = False 
            else:
                pass
        self._build_selection_buffer ()
        print ('selected atoms: ',len(self.selected_atoms))
    

    def selecting_by_molecule (self, selected_atom, disable = True):
        """       
        when disable = True
        If the object selection is disabled, all atoms within the 
        residue will be set to False 
        
        The disable variable does not allow, if the selected 
        atom is already in the selected list, to be removed. 
        
        The disable variable is "False" for when we use 
        selection by area (selection box)  """
        
        self._clear_selection_buffer(selected_atom)

            
        #if selected
        
        # if the selected atoms is not in the selected list
        #if selected_atom not in self.selected_atoms:
        #    
        #    for atom in selected_atom.residue.atoms:
        #        '''print (len(selected.residue.atoms), atom.name, atom.index)'''
        #        
        #        # the atom is not on the list -  add atom by atom
        #        if atom not in self.selected_atoms:
        #            self.selected_atoms.append(atom)
        #            atom.selected = True
        #        # the atom IS on the list - do nothing 
        #        else:
        #            pass
        #
        ## else: if the selected atoms IS in the selected list
        #else:
        #    if disable:
        #        # So, add all atoms  - selected residue <- selected.resi
        #        for atom in selected_atom.residue.atoms:
        #            
        #            # the atom is not on the list -  add atom by atom
        #            if atom in self.selected_atoms:
        #                index = self.selected_atoms.index(atom)
        #                self.selected_atoms.pop(index)
        #                atom.selected = False
        #
        #            # the atom IS in the list - do nothing 
        #            else:
        #                pass   
        #
        #    else:
        #        pass
        #   


    def unselecting_by_indexes (self, vismol_object = None, indexes = []):
        """ Function doc """
        print (indexes)
        #for atom in vismol_object.atoms: self.vismolSession.vismol_objects_dic.items()
        
        if vismol_object:
            for i in indexes:
                vismol_object.atoms[i].selected = False
        
        else:
            #for vismol_object in self.vismolSession.vismol_objects:
            for vobj_index, vismol_object in self.vismolSession.vismol_objects_dic.items():
                if indexes:
                    for i in indexes:
                        vismol_object.atoms[i].selected = False
                else:
                    for i in range(0, len(vismol_object.atoms)):
                        vismol_object.atoms[i].selected = False
                
                
        self._build_selection_buffer()
        #self.active = True
        self.build_selected_atoms_coords_and_selected_objects_from_selected_atoms ()

    def selecting_by_indexes (self, vismol_object = None, indexes = [], clear = False):
        """ Function doc """
        print (indexes)
        
        if clear:
            self._clear_selection_buffer ()
        #for atom in vismol_object.atoms:
        for i in indexes:
            vismol_object.atoms[i].selected = True
        
        self._build_selection_buffer()
        #self.active = True
        self.build_selected_atoms_coords_and_selected_objects_from_selected_atoms ()
    
    def invert_selection (self, vismol_object = None):
        """ not workign """
        if vismol_object:
            for atom in vismol_object.atoms:
                if atom.selected:
                    atom.selected = False
                else:
                    atom.selected = True
        
        else:
            #for vismol_object in self.vismolSession.vismol_objects:
            for index, vismol_object in self.vismolSession.vismol_objects_dic.items():
                for atom in vismol_object.atoms:
                    if atom.selected:
                        print (atom.name, atom.selected)
                        atom.selected = False
                    else:
                        atom.selected = True
            
        self._build_selection_buffer()
        #self.active = True
        self.build_selected_atoms_coords_and_selected_objects_from_selected_atoms ()
    
    def _build_selection_buffer (self):
        """ Function doc """
        self.selected_atoms = []
        
        #for Vobject in self.vismolSession.vismol_objects:
        for index, Vobject in self.vismolSession.vismol_objects_dic.items():
            for atom in Vobject.atoms:
                if atom.selected:
                    #print ('',atom )
                    self.selected_atoms.append(atom)
                else:
                    pass
        
        
    def _clear_selection_buffer (self,  selected_atom = None):
        ''' If the object selection is disabled, 
        all atoms in the system will be set to False '''

        #------------------------------------------------
        # Clearing the selections buffer
        # assigning false to all atoms loaded into memory
        if self.active:
            pass
        else:
            #for Vobject in self.vismolSession.vismol_objects:
            for index, Vobject in self.vismolSession.vismol_objects_dic.items():
                for atom in Vobject.atoms:
                    atom.selected = False
        #------------------------------------------------""

    def selection_function_viewing (self, selected, _type = None, disable =  True):
     
        '''
        Takes a selected atom and passes it to the appropriate selection function.
        '''        
        
        if _type:
            selection_mode2 = _type
        else:
            selection_mode2 = self._selection_mode
        
        if selected is None:
            self.selected_atoms = []
            self.selected_residues  = []
            self.active = False
        
        else:
            if selection_mode2 == 'atom':
                self.selecting_by_atom (selected, disable)
            
            elif selection_mode2 == 'residue':
                self.selecting_by_residue (selected, disable)

            elif selection_mode2 == 'chain':
                self.selecting_by_chain (selected, disable)
            else:
                pass
            
            self.active = True
        
        

        self.build_selected_atoms_coords_and_selected_objects_from_selected_atoms ()
    
    def build_selected_atoms_coords_and_selected_objects_from_selected_atoms (self):
        """ Function doc """
        

        self.selected_atoms_coords = []
        self.selected_objects          = {}
        for atom in self.selected_atoms:
            
            if atom.Vobject in self.selected_objects:
                self.selected_objects[atom.Vobject] += [atom.index-1]
            else:
                self.selected_objects[atom.Vobject] = [atom.index-1]
            
            #coords =  atom.coords()
            #self.selected_atoms_coords = self.selected_atoms_coords + coords
        
        
        for vobject in self.selected_objects:
            self.selected_objects[vobject] =  np.array(self.selected_objects[vobject], dtype=np.uint32)         
        
        #print  (self.selected_atoms_coords)      
        #for vobject in     self.selected_objects:
            ##print(vobject.name,self.selected_objects[vobject], 'selection_function_viewing' )
            #print('Selection defiend with  ',len(self.selected_atoms), 'atom(s)')
        
            
