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


class Config:
    """
    The Config class 
    """
    
    def __init__ (self):
        """ Class initialiser """
        pass
    def set_pk_sphr_selection_color (self, color, pk = 1):
        """ Function doc """
        
        if pk == 1:
            self.picking_selections.pk_scolor["pk1"] = color
        elif pk == 2:
            self.picking_selections.pk_scolor["pk2"] = color
        elif pk == 3:
            self.picking_selections.pk_scolor["pk3"] = color
        elif pk == 4:
            self.picking_selections.pk_scolor["pk4"] = color
        else:
            return False
        
        
    
    
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
    
    def __init__(self, vismol_widget: "VismolGTKWidget"=None, vismol_config: "VismolConfig"=None):
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
        # if vm_config:
        #     self.vm_config = vm_config
        # else:
        #     self.vm_config = VismolConfig()
        
        # self.toolkit = toolkit
        # if toolkit == "Gtk_3.0":
        #     from vismol.gui.vismol_widget_main import VismolWidgetMain
        #     # TODO: The following three lines should be improved. Right now
        #     #       there is a cyclic dependency between VismolWidgetMain,
        #     #       VismolGLMain and VismolSession. This approach is working
        #     #       now, but it may cause problems in future versions :(
        #     self.vm_widget = VismolWidgetMain(self.vm_config)
        #     self.vm_widget.vm_session = self
        #     self.vm_widget.vm_glcore.vm_session = self
        
        #     self.selection_box_frame = None
        #     self.vm_glcore = self.vm_widget.vm_glcore
        #     self.vm_glcore.queue_draw()
        # elif toolkit == "Qt5":
        #     self.vm_widget = None
        #     logger.error("Not implemented yet for Qt5 :(")
        #     raise NotImplementedError("Not implemented yet for Qt5 :(")
        #     exit(1)
        # else:
        #     self.vm_widget = None
        #     logger.error("Toolkit not defined or syntax error, try 'Gtk_3.0'. Quitting.")
        #     raise RuntimeError("Toolkit not defined or syntax error, try 'Gtk_3.0'. Quitting.")
        #     exit(1)
        
        # TODO: Add support for Qt5? This parameter is only used for the
        #       center_on_coordinates function in VismolGLMain class.
        # self.toolkit = "Gtk_3.0"
        # self.vm_config = vismol_config
        # self.vm_widget = vismol_widget
        # self.selection_box_frame = None
        # self.vm_glcore = self.vm_widget.vm_glcore
        # self.vm_glcore.queue_draw()
        # self.frame = 0
        # self.vm_objects_dic = {}
        # self.atom_dic_id = {}
        # self.atom_id_counter = np.uint32(0)
        self.picking_selection_mode = False
        self.selections = {"sel_00": VMSele(self)}
        self.current_selection = "sel_00"
        self.picking_selections = VMPick(self)
        # self.periodic_table = PeriodicTable()
        self.vm_geometric_object_dic = {"pk1pk2":None, "pk2pk3":None,
                "pk3pk4":None, "pk1":None, "pk2":None, "pk3":None, "pk4":None}
    
    def _add_vismol_object(self, vm_object, show_molecule=True, autocenter=True):
        """ Function doc """
        
        # Check if the vismol_object's index already exists in the dictionary
        if vm_object.index in self.vm_objects_dic.keys():
            logger.warning("The VismolObject with id {} already exists. \
                The data will be overwritten.".format(vm_object.index))
        
        # Add the vismol_object to the dictionary with a new index
        self.vm_objects_dic[len(self.vm_objects_dic)] = vm_object
        
        # Update the atom_id_counter by adding the number of atoms in the new vismol_object
        # TODO: Is this really necessary? The number of atoms can be obtained
        #       from the atom_dictionary 'self.atom_dic_id'
        self.atom_id_counter += len(vm_object.molecule.atoms)
        
        # Iterate through the atoms in the vismol_object and add them to the atom_dic_id dictionary
        # using their unique_id as the key and the atom object as the value.
        for atom in vm_object.molecule.atoms.values():
            self.atom_dic_id[atom.unique_id] = atom
        
        # If show_molecule is True, create representations for "lines" and "nonbonded" to display the object.
        if show_molecule:
            vm_object.create_representation(rep_type="lines")
            vm_object.create_representation(rep_type="nonbonded")
            # TODO: Add support for cartoon representation
            # vm_object.create_representation(rep_type="cartoon")
            # TODO: Add support for surfae representation
            # vm_object.create_representation(rep_type="surface")
            
            # If autocenter is True, center the view on the mass center of the vismol_object.
            if autocenter:
                self.vm_glcore.center_on_coordinates(vm_object, vm_object.molecule.geom_center)
            else:
                self.vm_glcore.queue_draw()
    
    def load_molecule(self, infile):
        """ Probably would be better to join this with _add_vismol_object
        """
        vm_object, show_molecule = parser.parse_file(self, infile)
        vm_object.set_model_matrix(self.vm_glcore.model_mat)
        vm_object.generate_color_vectors()
        # TODO: Check where this is used vm_object.active
        vm_object.active = True
        self._add_vismol_object(vm_object, show_molecule=show_molecule)
    
    # def _change_attributes_for_atoms(self, atoms, rep_type, show):
    #     """ Function doc """
    #     for atom in atoms:
    #         try:
    #             if show:
    #                 setattr(atom, rep_type, True)
    #             else:
    #                 setattr(atom, rep_type, False)
    #         except AttributeError as ae:
    #             logger.error("Representation of type {} not implemented".format(rep_type))
    #             logger.error(ae)
    
    def show_or_hide(self, rep_type="lines", selection=None, show=True):
        """ Function doc """
        if selection is None:
            selection = self.selections[self.current_selection]
        
        self._change_attributes_for_atoms(selection.selected_atoms, rep_type, show)
        for vm_object in selection.selected_objects:
            if vm_object.representations[rep_type] is None:
                vm_object.create_representation(rep_type = rep_type)
            show_hide_indexes = []
            if rep_type == "lines":
                for bond in vm_object.molecule.bonds:
                    if bond.atom_i.lines and bond.atom_j.lines:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "sticks":
                for bond in vm_object.molecule.bonds:
                    if bond.atom_i.sticks and bond.atom_j.sticks:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dash":
                for bond in vm_object.molecule.bonds:
                    if bond.atom_i.dash and bond.atom_j.dash:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dynamic":
                self.define_dynamic_bonds(selection = selection)
                for bond in vm_object.molecule.bonds:
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
                for bond in vm_object.molecule.bonds:
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
                #print('show_hide_indexes',show_hide_indexes)
                #logger.error("Not implementer for 'ribbon' yet.")
                #raise NotImplementedError("Not implementer for 'ribbon' yet.")
            
            
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
                print('vm_object.representations[rep_type].active = True')
                vm_object.representations[rep_type].define_new_indexes_to_vbo(show_hide_indexes)
                vm_object.representations[rep_type].active = True
                vm_object.representations[rep_type].was_rep_ind_modified = True
                vm_object.representations[rep_type].was_sel_ind_modified = True
            else:
                vm_object.representations[rep_type].active = False
                print('vm_object.representations[rep_type].active = False')
        
        self.vm_widget.queue_draw()
        return selection
    
    def forward_frame(self):
        """ Function doc """
        #return 0
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
        
        #if self.picking_selection_mode:
        #    self.picking_selections.print_pk_distances()
    
    def reverse_frame(self):
        """ Function doc """
        if self.frame - 1 >= 0:
            self.frame -= 1
            self.vm_glcore.updated_coords = True
        else:
            self.vm_glcore.updated_coords = False
        
        #if self.picking_selection_mode:
        #    self.picking_selections.print_pk_distances()
    
    def set_frame(self, frame=0):
        """ Function doc """
        assert frame >= 0
        self.frame = np.uint32(frame)
        
        #self.picking_selections
        #for 
        self.picking_selections.update_pki_pkj_rep_coordinates()
        self.vm_widget.queue_draw()
        print('\n\n\nuhuuu')
        #if self.picking_selection_mode:
        #    self.picking_selections.print_pk_distances()
    
    def get_frame(self):
        """ Function doc """
        return self.frame
    
    def _selection_function_set(self, selected, _type=None, disable=True):
        """ Function doc """
        if self.picking_selection_mode: # True for picking mode
            if selected:
                assert len(selected) == 1
                selected = list(selected)[0]
                self.picking_selections.selection_function_picking(selected)
            else:
                self.picking_selections.selection_function_picking(None)
        else: # False for viewing mode
            self.selections[self.current_selection].selection_function_viewing_set(selected, _type, disable)
    
    def viewing_selection_mode(self, sel_type="atom"):
        """ Function doc """
        if self.selection_box_frame:
            self.selection_box_frame.change_sel_type_in_combobox(sel_type)
        self.selections[self.current_selection].selection_mode = sel_type
    
    def define_dynamic_bonds(self, selection = None):
        """ Function doc """
        if selection:
            pass
        else:
            selection = self.selections[self.current_selection]
        
        selection_dict = {}
        vobject = None
        for atom in selection.selected_atoms:
            selection_dict[atom.atom_id] = atom
            vobject = atom.vm_object
        
        tolerance = self.vm_config.gl_parameters['bond_tolerance']
        vobject.dynamic_bonds = []
        for frame in range(len(vobject.frames)):            
            bonds = vobject.build_bonded_and_nonbonded_atoms(selection=selection_dict, frame=frame, tolerance = tolerance)
            vobject.dynamic_bonds.append(bonds)
            #print(len(bonds), bonds)
        #print(vobject.dynamic_bonds)

    # def hide_axes (self):
    #     """ Function doc """
    #     self.vm_glcore.show_axis = False
    # def show_axes (self):
    #     """ Function doc """
    #     self.vm_glcore.show_axis = True

        
    
    # def delete_vismol_object_by_index(self, index):
    #     """ Function doc
    #     """
    #     try:
    #         vm_object = self.vm_objects_dic.pop(index)
    #     except KeyError:
    #         logger.warning("VismolObject with index {} not found".format(index))
    #         return None
    #     return vm_object
    
    # def select(self, vismol_object=None, indexes=None, sele=None):
    #     """ Function doc """
    #     if vismol_object is None:
    #         vismol_object = self.vm_objects[-1]
        
    #     if sele is None:
    #         sele = self.current_selection
        
    #     if indexes == "all":
    #         self.selections[sele].selecting_by_indexes(vismol_object=vismol_object,
    #                                                    indexes=range(0, int(len(vismol_object.atoms)/2)))
    #     self.vm_widget.queue_draw()
    
    # def disable_by_index(self, index):
    #     """ When the variable "dictionary" is active, the function accesses 
    #         a vismol object through the dictionary "self.vm_objects_dic". 
    #         Each vismol object has a unique access key (int), which, in 
    #         easyhybrid, is generated in the method: add_vismol_object.
            
    #         In the vismol interface the enable_by_index/disable_by_index methods
    #         access the vismol objects by their position in the "self.vm_objects" 
    #         list (this is because when an object is deleted in the vismol 
    #         interface, the treeview"s liststore is rewritten)
    #     """
    #     try:
    #         self.vm_objects_dic[index].active = False
    #         self.vm_glcore.queue_draw()
    #     except KeyError:
    #         logger.error("VismolObject with index {} not found".format(index))
    #         return False
    #     return True
    
    # def enable_by_index(self, index):
    #     """ When the variable "dictionary" is active, the function accesses 
    #         a vismol object through the dictionary "self.vm_objects_dic". 
    #         Each vismol object has a unique access key (int), which, in 
    #         easyhybrid, is generated in the method: add_vismol_object.
            
    #         In the vismol interface the enable_by_index/disable_by_index methods
    #         access the vismol objects by their position in the "self.vm_objects" 
    #         list (this is because when an object is deleted in the vismol 
    #         interface, the treeview"s liststore is rewritten)
    #     """
    #     try:
    #         self.vm_objects_dic[index].active = True
    #         self.vm_glcore.queue_draw()
    #     except KeyError:
    #         logger.error("VismolObject with index {} not found".format(index))
    #         return False
    #     return True
    
    # def edit_by_index(self, index):
    #     """ Function doc
    #     """
    #     try:
    #         self.vm_objects_dic[index].moving = not self.vm_objects_dic[index].moving
    #         self.vm_glcore.queue_draw()
    #     except KeyError:
    #         logger.error("VismolObject with index {} not found".format(index))
    #         return False
    #     return True
    
    # def set_color_by_index(self, vismol_object, indexes=None, color=None):
    #     """ NOT SURE WHAT THIS FUNCTION DOES
    #     """
    #     if indexes is None:
    #         indexes = []
    #     if color is None:
    #         color = np.array([0.9, 0.9, 0.9], dtype=np.float32)
        
    #     for atom_index in indexes:
    #         vismol_object.atoms[atom_index].color = color
    #     vismol_object.generate_color_vectors(do_colors=True, do_colors_idx=False,
    #                                           do_colors_raindow=False, do_vdw_dot_sizes=False,
    #                                           do_cov_dot_sizes=False)
    #     self.vm_widget.queue_draw()
    #     for rep  in vismol_object.representations.keys():
    #         if vismol_object.representations[rep]:
    #             vismol_object.representations[rep]._set_colors_to_buffer()
    #             # try:
    #             #     vismol_object.representations[rep]._set_colors_to_buffer()
    #             # except:
    #             #     print('"VisMol/vModel/Representations.py, line 123, in _set_colors_to_buffer GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ctypes.ArgumentError: argument 2: <class "TypeError">: wrong type"')
    #     return True
    
