#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import numpy as np
from logging import getLogger
from vismol.utils import parser
from vismol.model.atom import Atom
from vismol.core.vismol_object import VismolObject
from vismol.core.vismol_session import VismolSession


logger = getLogger(__name__)


class VismolSessionBuilder(VismolSession):
    """
    The VismolSessionMain class acts as a bridge between Vismol objects containing
    the molecular structures, the OpenGL framework and the GTK library.
    
    Attributes:
        atom_dic_id (Dict): Dictionary containing all the atoms in the session,
                regardless of the VismolObject that they belong to.
        atom_id_counter (int): Integer representing the number of atoms present
                in the session, regardless of the VismolObject that they belong
                to. This value cannot be the same as len(atom_dic_id), since it
                may be cases in that some atoms could be removed from the
                session. TODO: Check if this value is utilized. TO_REVIEW
        current_selection (str): The current selection active. This corresponds
                to the key in the selections dictionary.
        frame (int): Te current frame that the session is in. It depends of the
                VismolObject elements in the session, it should be set to the
                VismolObject selected. TODO: Implement similar to VMD?
        periodic_table (PeriodicTable): The PeriodicTable object contains the
                information for the colors of the elements and also some of
                their atomic properties.
        picking_selection_mode (bool): Whether the click on the window is to
                perform movement actions or picking actions, if True, then picking.
        picking_selections (VismolPickingSelection): The VismolPickingSelection
                object for this session. TODO: Add description here. TO_REVIEW
        selections (Dict): Dictionary of VismolViewingSelection objects.
                TODO: Add description here. TO_REVIEW
        toolkit (str): The window manager library employed to build the GUI.
                Currently the only supported library is Gtk 3.0. TODO: Add
                support for Qt5 and Tkinter?
        vm_config (VismolConfig): The configuration object containing all the
                parameters used for the rendering and layout of the window.
        vm_geometric_object_dic (Dict): Dictionary containing the atoms selected
                with the picking mode. This atoms are used to calculate some
                properties like distance and angles.
        vm_object_dic (Dict): Dictionary containing all the VismolObject elements
                in the session.
        widget (object): Depending on the toolkit selected, this object should
                implement all the functions that connect the Vismol functionality
                with the user interface. For the Gtk library, it corresponds to
                a VismolWidgetMain object.
    """
    
    def __init__(self, vismol_widget: "VismolGTKWidget", vismol_config: "VismolConfig"):
        """
        Args:
            toolkit (str): The window manager library employed to build the GUI.
                    Currently the only supported library is Gtk 3.0. TODO: Add
                    support for Qt5 and Tkinter?
            widget (object): Depending on the toolkit selected, this object
                    should implement all the functions that connect the Vismol
                    functionality with the user interface. For the Gtk library,
                    it corresponds toa VismolWidgetMain object.
            vm_config (VismolConfig): The configuration object containing all the
                    parameters used for the rendering and layout of the window.
        """
        
        super(VismolSessionBuilder, self).__init__(vismol_widget, vismol_config)
        """
        """
        self.atom_ids = 0
    
    def add_new_atom(self, coords: np.array) -> None:
        if len(self.vm_objects_dic) == 0:
            vm_object = VismolObject(self, index=0, name="Builder")
            vm_object.set_model_matrix(self.vm_widget.vm_glcore.model_mat)
            atom = Atom(vm_object, name="C", index=self.atom_ids,
                        atom_id=self.atom_ids, coords=coords)
            atom.spheres = True
            atom.sticks = True
            atom.unique_id = self.atom_ids
            self.atom_dic_id[atom.unique_id] = atom
            vm_object.add_new_atom(atom, self.atom_ids)
            vm_object.generate_color_vectors()
            vm_object.active = True
            self.vm_objects_dic[0] = vm_object
        else:
            atom = Atom(self.vm_objects_dic[0], name="C", index=self.atom_ids,
                        atom_id=self.atom_ids, coords=coords)
            atom.spheres = True
            atom.sticks = True
            atom.unique_id = self.atom_ids
            self.atom_dic_id[atom.unique_id] = atom
            self.vm_objects_dic[0].add_new_atom(atom, self.atom_ids)
            self.vm_objects_dic[0].generate_color_vectors()
            self.vm_objects_dic[0].active = True
        
        self.atom_ids += 1
        self.vm_objects_dic[0].create_representation(rep_type="spheres")
        self.vm_objects_dic[0].create_representation(rep_type="sticks")
        self.vm_objects_dic[0].representations["spheres"].define_new_indexes_to_vbo(list(self.atom_dic_id.keys()))
        self.vm_objects_dic[0].representations["sticks"].define_new_indexes_to_vbo(list(self.atom_dic_id.keys()))
        self.vm_objects_dic[0].representations["spheres"].was_rep_ind_modified = True
        self.vm_objects_dic[0].representations["sticks"].was_rep_ind_modified = True
        
        self.vm_glcore.queue_draw()
    

