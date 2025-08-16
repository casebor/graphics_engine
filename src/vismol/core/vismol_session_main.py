#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import numpy as np
from logging import getLogger
from vismol.utils import parser
from vismol.utils.elements import PeriodicTable
from vismol.core.vismol_session import VismolSession
from vismol.core.vismol_selections import VismolPickingSelection as VMPick
from vismol.core.vismol_selections import VismolViewingSelection as VMSele


logger = getLogger(__name__)


class VismolSessionMain(VismolSession):
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
        
        super(VismolSessionMain, self).__init__(vismol_widget, vismol_config)
        self.picking_selection_mode = False
        self.selections = dict({"sel_00": VMSele(self)})
        self.current_selection = "sel_00"
        self.picking_selections = VMPick(self)
        self.vm_geometric_object_dic = {"pk1":None, "pk2":None, "pk3":None,
                    "pk4":None, "pk1pk2":None, "pk2pk3":None, "pk3pk4":None}
    
    def _add_vismol_object(self, vismol_object: "VismolObject",
                           show_molecule: bool=True,
                           autocenter: bool=True) -> None:
        """ Function doc """
        
        # Check if the vismol_object's index already exists in the dictionary
        if vismol_object.index in self.vm_objects_dic.keys():
            logger.warning("The VismolObject with id {} already exists. \
                The data will be overwritten.".format(vismol_object.index))
        
        # Add the vismol_object to the dictionary with a new index
        self.vm_objects_dic[len(self.vm_objects_dic)] = vismol_object
        
        # Update the atom_id_counter by adding the number of atoms in the new vismol_object
        # TODO: Is this really necessary? The number of atoms can be obtained
        #       from the atom_dictionary 'self.atom_dic_id'
        self.atom_id_counter += len(vismol_object.molecule.atoms)
        
        # Iterate through the atoms in the vismol_object and add them to the atom_dic_id dictionary
        # using their unique_id as the key and the atom object as the value.
        for atom in vismol_object.molecule.atoms.values():
            self.atom_dic_id[atom.unique_id] = atom
        
        # If show_molecule is True, create representations for "lines" and "nonbonded" to display the object.
        if show_molecule:
            vismol_object.create_representation(rep_type="lines")
            vismol_object.create_representation(rep_type="nonbonded")
            # TODO: Add support for cartoon representation
            # vismol_object.create_representation(rep_type="cartoon")
            # TODO: Add support for surfae representation
            # vismol_object.create_representation(rep_type="surface")
            
            # If autocenter is True, center the view on the mass center of the vismol_object.
            if autocenter:
                self.vm_glcore.center_on_coordinates(vismol_object, vismol_object.molecule.get_geometry_center())
            else:
                self.vm_glcore.queue_draw()
    
    def _selection_function_set(self, selected: set, disable: bool=True) -> None:
        """ Function doc """
        # TODO: the selected parameter should be an integer, since it is
        #       only one ID. Discuss it with Bachega
        if self.picking_selection_mode: # True for picking mode
            if selected:
                assert len(selected) == 1
                selected = list(selected)[0]
                self.picking_selections.selection_function_picking(selected)
            else:
                self.picking_selections.selection_function_picking(None)
        else: # False for viewing mode
            self.selections[self.current_selection].selection_function_viewing_set(selected, disable)
    
    def define_dynamic_bonds(self, vmv_selection: VMSele=None) -> None:
        """ Function doc """
        if vmv_selection:
            pass
        else:
            vmv_selection = self.selections[self.current_selection]
        
        selection_dict = {}
        vm_object = None
        for atom in vmv_selection.selected_atoms:
            selection_dict[atom.atom_id] = atom
            vm_object = atom.vm_object
        
        tolerance = self.vm_config.gl_parameters['bond_tolerance']
        vm_object.dynamic_bonds = []
        for frame in range(len(vm_object.frames)):            
            bonds = vm_object.build_bonded_and_nonbonded_atoms(selection=selection_dict, frame=frame, tolerance = tolerance)
            vm_object.dynamic_bonds.append(bonds)
    
    def forward_frame(self) -> None:
        """ Function doc """
        frame = self.frame + 1
        for i, vm_object in enumerate(self.vm_objects_dic.values()):
            if frame < vm_object.molecule.frames.shape[0]:
                self.frame += 1
                self.vm_glcore.updated_coords = True
                break
            else:
                pass
        else:
            self.vm_glcore.updated_coords = False
    
    def get_frame(self) -> int:
        """ Function doc """
        return self.frame
    
    def load_molecule(self, infile: str) -> None:
        """ Probably would be better to join this with _add_vismol_object
        """
        vm_object, show_molecule = parser.parse_file(self, infile)
        vm_object.set_model_matrix(self.vm_glcore.model_mat)
        vm_object.generate_color_vectors()
        # TODO: Check where this is used vm_object.active
        vm_object.active = True
        self._add_vismol_object(vm_object, show_molecule=show_molecule)
    
    def reverse_frame(self) -> None:
        """ Function doc """
        if self.frame - 1 >= 0:
            self.frame -= 1
            self.vm_glcore.updated_coords = True
        else:
            self.vm_glcore.updated_coords = False
    
    def set_frame(self, frame: int=0) -> None:
        """ Function doc """
        assert frame >= 0
        self.frame = np.uint32(frame)
        self.picking_selections.update_pki_pkj_rep_coordinates()
        self.vm_widget.queue_draw()
    
    def show_or_hide(self, rep_type: str="lines", vmv_selection: VMSele=None,
                     show: bool=True) -> VMSele:
        """ Function doc """
        if vmv_selection is None:
            vmv_selection = self.selections[self.current_selection]
        
        self._change_attributes_for_atoms(vmv_selection.selected_atoms, rep_type, show)
        for vm_object in vmv_selection.selected_objects:
            if vm_object.representations[rep_type] is None:
                vm_object.create_representation(rep_type = rep_type)
            show_hide_indexes = []
            if rep_type == "lines":
                for bond in vm_object.molecule.topology.bonds:
                    if bond.atom_i.lines and bond.atom_j.lines:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "sticks":
                for bond in vm_object.molecule.topology.bonds:
                    if bond.atom_i.sticks and bond.atom_j.sticks:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dash":
                for bond in vm_object.molecule.topology.bonds:
                    if bond.atom_i.dash and bond.atom_j.dash:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dynamic":
                self.define_dynamic_bonds(selection = vmv_selection)
                for bond in vm_object.molecule.topology.bonds:
                    if bond.atom_i.dynamic and bond.atom_j.dynamic:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dots":
                for atom in vm_object.molecule.atoms.values():
                    if atom.dots:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "nonbonded":
                for atom in vm_object.molecule.atoms.values():
                    if atom.nonbonded:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "impostor":
                for atom in vm_object.molecule.atoms.values():
                    if atom.impostor:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "spheres":
                for atom in vm_object.molecule.atoms.values():
                    if atom.spheres:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "stick_spheres":
                for bond in vm_object.molecule.topology.bonds:
                    if bond.atom_i.stick_spheres and bond.atom_j.stick_spheres:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
                
                
                #for atom in vm_object.molecule.atoms.values():
                #    if atom.stick_spheres:
                #        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "vdw_spheres":
                for atom in vm_object.molecule.atoms.values():
                    if atom.vdw_spheres:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "ribbons": # add by bachega at 02/02/2023
                for bond in vm_object.molecule.c_alpha_bonds:
                    if bond.atom_i.ribbons and bond.atom_j.ribbons:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "ribbon_sphere":
                for atom in vm_object.molecule.atoms.values():
                    if atom.ribbon_sphere and atom.name=='CA' and atom.residue.is_protein:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "surface":
                logger.error("Not implementer for 'surface' yet.")
                raise NotImplementedError("Not implementer for 'surface' yet.")
            
            elif rep_type == "cartoon":
                logger.error("Not implementer for 'cartoon' yet.")
                raise NotImplementedError("Not implementer for 'cartoon' yet.")
            
            elif rep_type == "labels":
                for atom in vm_object.molecule.atoms.values():
                    if atom.labels:
                        show_hide_indexes.append(atom.atom_id)
            
            if len(show_hide_indexes) > 0:
                # print('vm_object.representations[rep_type].active = True')
                vm_object.representations[rep_type].define_new_indexes_to_vbo(show_hide_indexes)
                vm_object.representations[rep_type].active = True
                vm_object.representations[rep_type].was_rep_ind_modified = True
                vm_object.representations[rep_type].was_sel_ind_modified = True
            else:
                vm_object.representations[rep_type].active = False
                # print('vm_object.representations[rep_type].active = False')
        
        self.vm_widget.queue_draw()
        return vmv_selection
    
    def viewing_selection_mode(self, sel_type: str="atom") -> None:
        """ Function doc """
        if self.selection_box_frame:
            self.selection_box_frame.change_sel_type_in_combobox(sel_type)
        self.selections[self.current_selection].selection_mode = sel_type
    