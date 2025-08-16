#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import numpy as np
# from typing import List
from logging import getLogger
from vismol.utils.elements import PeriodicTable


logger = getLogger(__name__)


class VismolSession:
    """
    The VismolSession class acts as a bridge between Vismol objects containing
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
        
        # TODO: Add support for Qt5? This parameter is only used for the
        #       center_on_coordinates function in VismolGLMain class.
        self.toolkit = "Gtk_3.0"
        self.vm_config = vismol_config
        self.vm_widget = vismol_widget
        self.selection_box_frame = None
        self.vm_glcore = self.vm_widget.vm_glcore
        self.vm_glcore.queue_draw()
        self.frame = 0
        self.vm_objects_dic = {}
        self.atom_dic_id = {}
        # TODO: Is this variable needed? This number can be obtained from the
        #       len(self.vm_objects_dic)
        self.atom_id_counter = np.uint32(0)
        # self.picking_selection_mode = False
        # self.selections = {"sel_00": VMSele(self)}
        # self.current_selection = "sel_00"
        # self.picking_selections = VMPick(self)
        self.periodic_table = PeriodicTable()
        # self.vm_geometric_object_dic = {"pk1pk2":None, "pk2pk3":None,
        #         "pk3pk4":None, "pk1":None, "pk2":None, "pk3":None, "pk4":None}
    
    def _change_attributes_for_atoms(self, atoms: list, rep_type: str,
                                     show: bool) -> None:
        """ Function doc """
        for atom in atoms:
            try:
                if show:
                    setattr(atom, rep_type, True)
                else:
                    setattr(atom, rep_type, False)
            except AttributeError as ae:
                logger.error("Representation type {} not implemented".format(rep_type))
                logger.error(ae)
    
    def hide_axis(self) -> None:
        """ Function doc """
        self.vm_glcore.show_axis = False
    
    def show_axis(self) -> None:
        """ Function doc """
        self.vm_glcore.show_axis = True
    
    def _add_vismol_object(self, vismol_object: "VismolObject",
                           show_molecule: bool, autocenter: bool) -> None:
        logger.critical("NotImplementedError, the child class must implement _add_vismol_object")
        raise NotImplementedError("Subclasses must implement this method")
    
    def load_molecule(self, infile: str) -> None:
        logger.critical("NotImplementedError, the child class must implement load_molecule")
        raise NotImplementedError("Subclasses must implement this method")
    
    def show_or_hide(self, rep_type: str, vmv_selection: "VismolViewingSelection",
                     show: bool) -> "VismolViewingSelection":
        logger.critical("NotImplementedError, the child class must implement show_or_hide")
        raise NotImplementedError("Subclasses must implement this method")
    
    def forward_frame(self) -> None:
        logger.critical("NotImplementedError, the child class must implement forward_frame")
        raise NotImplementedError("Subclasses must implement this method")
    
    def reverse_frame(self) -> None:
        logger.critical("NotImplementedError, the child class must implement reverse_frame")
        raise NotImplementedError("Subclasses must implement this method")
    
    def set_frame(self, frame: int=0) -> None:
        logger.critical("NotImplementedError, the child class must implement set_frame")
        raise NotImplementedError("Subclasses must implement this method")
    
    def get_frame(self) -> int:
        logger.critical("NotImplementedError, the child class must implement get_frame")
        raise NotImplementedError("Subclasses must implement this method")
    
    def _selection_function_set(self, selected: set, disable: bool) -> None:
        # TODO: the selected parameter should be an integer, since it is
        #       only one ID. Discuss it with Bachega
        logger.critical("NotImplementedError, the child class must implement _selection_function_set")
        raise NotImplementedError("Subclasses must implement this method")
    
    def viewing_selection_mode(self, sel_type: str) -> None:
        logger.critical("NotImplementedError, the child class must implement viewing_selection_mode")
        raise NotImplementedError("Subclasses must implement this method")
    
