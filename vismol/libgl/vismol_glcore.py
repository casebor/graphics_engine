#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_glcore.py
#  
#  Copyright 2022 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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

import time
import numpy as np
from OpenGL import GL
from libgl.glaxis import GLAxis
from libgl.glcamera import GLCamera
from libgl.vismol_font import VismolFont
from libgl.dynamic_line import DynamicLine
from libgl.selection_box import SelectionBox
import libgl.shapes as shapes
import libgl.shaders.pick as shaders_pick
import libgl.shaders.dots as shaders_dots
import libgl.shaders.lines as shaders_lines
import libgl.shaders.wires as shaders_wires
import libgl.shaders.sticks as shaders_sticks
import libgl.shaders.cartoon as shaders_cartoon
import libgl.shaders.surface as shaders_surface
import libgl.shaders.spheres as shaders_spheres
import libgl.shaders.impostor as shaders_impostor
import libgl.shaders.nonbonded as shaders_nonbonded
import libgl.shaders.vm_freetype as shaders_vm_freetype
import utils.matrix_operations as mop


class VismolGLCore:
    
    def __init__(self, widget, vismol_session=None, width=640.0, height=420.0):
        """ Constructor of the class.
            
            Keyword arguments:
            vismol_session - 
        """
        self.parent_widget = widget
        self.vm_session = vismol_session
        self.vm_config = self.vm_session.vm_config
        self.width = np.float32(width)
        self.height = np.float32(height)
        self.shader_programs = {}
        self.representations_available = {"dots", "lines", "nonbonded", "dotted_lines",
                                          "ribbon", "sticks", "spheres", "impostor",
                                          "surface", "cartoon", "freetype",
                                          "picking_dots"}
    
    def initialize(self):
        """ Enables the buffers and other charasteristics of the OpenGL context.
            sets the initial projection, view and model matrices
            
            self.flag -- Needed to only create one OpenGL program, otherwise a bunch of
                         programs will be created and use system resources. If the OpenGL
                         program will be changed change this value to True
        """
        self.model_mat = np.identity(4, dtype=np.float32) # Not sure if this is used :S
        self.zero_reference_point = np.zeros(3, dtype=np.float32)
        self.glcamera = GLCamera(self.vm_config.gl_parameters["field_of_view"],
                                 self.width / self.height,
                                 np.array([0,0,10], dtype=np.float32),
                                 self.zero_reference_point)
        self.vm_font = VismolFont(color=[1, 1, 1, 0.6])
        self.axis = GLAxis()
        self.selection_box = SelectionBox()
        self.dynamic_line = DynamicLine()
        self.parent_widget.set_has_depth_buffer(True)
        self.parent_widget.set_has_alpha(True)
        self.frame = self.vm_session.frame
        self.scroll = self.vm_config.gl_parameters["scroll_step"]
        self.bckgrnd_color = self.vm_config.gl_parameters["background_color"]
        #                       Light Parameters                                
        self.light_position = self.vm_config.gl_parameters["light_position"]
        self.light_ambient_coef = self.vm_config.gl_parameters["light_ambient_coef"]
        self.light_shininess = self.vm_config.gl_parameters["light_shininess"]
        self.light_intensity = self.vm_config.gl_parameters["light_intensity"]
        #                              Variables                                
        self.right = self.width / self.height
        self.left = -self.right
        self.top = np.float32(1.0)
        self.bottom = np.float32(-1.0)
        self.button = None
        self.dist_cam_zrp = np.linalg.norm(self.glcamera.get_position())
        self.shader_flag = True
        self.modified_data = False
        self.modified_view = False
        self.dragging = False
        self.editing_mols = False
        self.show_axis = True
        self.ctrl = False
        self.shift = False
        self.atom_picked = None
        self.selection_box_picking = False
        self.picking = False
        self.picking_x, self.picking_y = None, None
        self.show_dynamic_line = False
        self.show_dynamic_line_x, self.show_dynamic_line_y = None, None
        self.show_selection_box = False
        self.show_selection_box_x, self.show_selection_box_y = None, None
        self.mouse_x, self.mouse_y = np.float32(0.0), np.float32(0.0)
        self.mouse_rotate, self.mouse_zoom, self.mouse_pan = False, False, False
        self.drag_pos_x, self.drag_pos_y, self.drag_pos_z = None, None, None
        self.chars, self.xyz_pos, self.uv_coords = None, None, None
    
    def resize_window(self, width, height):
        """ Resizing function, takes the widht and height of the widget
            and modifies the view in the camera acording to the new values
        
            Keyword arguments:
            width -- Actual width of the window
            height -- Actual height of the window
        """
        self.width = np.float32(width)
        self.height = np.float32(height)
        self.left = -self.width / self.height
        self.right = -self.left
        self.center_x = self.width / 2.0
        self.center_y = self.height / 2.0
        self.glcamera.viewport_aspect_ratio = self.width / self.height
        _proj_mat = mop.my_glPerspectivef(self.glcamera.field_of_view,
                                          self.glcamera.viewport_aspect_ratio,
                                          self.glcamera.z_near, self.glcamera.z_far)
        self.glcamera.set_projection_matrix(_proj_mat)
    
    def mouse_pressed(self, button_number, mouse_x, mouse_y):
        """ Function doc
        """
        left   = np.int32(button_number) == 1
        middle = np.int32(button_number) == 2
        right  = np.int32(button_number) == 3
        self.mouse_rotate = left   and not (middle or right)
        self.mouse_zoom   = right  and not (middle or left)
        self.mouse_pan    = middle and not (right  or left)
        self.mouse_x = np.float32(mouse_x)
        self.mouse_y = np.float32(mouse_y)
        self.drag_pos_x, self.drag_pos_y, self.drag_pos_z = self._mouse_pos(self.mouse_x, self.mouse_y)
        self.dragging = False
        if left:
            if self.shift:
                if self.show_dynamic_line:
                    self.dynamic_line.start = self.get_viewport_pos(mouse_x, mouse_y)
                    self.dynamic_line.end = self.get_viewport_pos(mouse_x, mouse_y)
                    self.dynamic_line.update_points()
                    self.dynamic_line_x = mouse_x
                    self.dynamic_line_y = self.height - mouse_y
                else:
                    self.show_selection_box = True
                    self.selection_box.start = self.get_viewport_pos(mouse_x, mouse_y)
                    self.selection_box.end = self.get_viewport_pos(mouse_x, mouse_y)
                    self.selection_box.update_points()
                    self.selection_box_x = mouse_x
                    self.selection_box_y = self.height - mouse_y
        if middle:
            self.picking_x = np.float32(mouse_x)
            self.picking_y = np.float32(mouse_y)
            self.picking = True
            self.queue_draw()
        if right:
            self.picking_x = np.float32(mouse_x)
            self.picking_y = np.float32(mouse_y)
            self.picking = True
            self.queue_draw()
    
    def mouse_released(self, button_number, mouse_x, mouse_y):
        """ Function doc
        int(event.button)
        
        info      = menu header info
        menu_type = "pick_menu" / "bg_menu" / "sele_menu" / "ob_menu"
        
        """
        left   = np.int32(button_number) == 1
        middle = np.int32(button_number) == 2
        right  = np.int32(button_number) == 3
        self.mouse_rotate = False
        self.mouse_zoom = False
        self.mouse_pan = False
        if self.dragging:
            if left:
                if self.shift:
                    self.selection_box_picking = True
                    self.show_selection_box = False
                    self.selection_box.start = None
                    self.selection_box.end = None
                    self.queue_draw()
        else:
            if left:
                self.picking_x = np.float32(mouse_x)
                self.picking_y = np.float32(mouse_y)
                self.picking = True
                self.button = 1
                self.queue_draw()
            if middle:
                if self.atom_picked is not None:
                    self.button = 2
                    self.center_on_atom(self.atom_picked)
                    self.atom_picked = None
            if right:
                # The right button (button = 3) always opens one of the available menus.
                self.button = 3
                menu_type = None
                # Check if there is anything in the selection list
                # If {} means that there are no selection points on the screen
                # Checks if vismol_session.current_selection has any selection.
                # Also needs to check whether "picking" mode is enabled.
                if not bool(self.vm_session.selections[self.vm_session.current_selection].selected_objects) \
                   or self.vm_session._picking_selection_mode:
                    # Checks if the list of atoms selected by the picking function has any elements. 
                    # If the list is empty, the pick menu is not shown.
                    if self.vm_session._picking_selection_mode \
                       and self.vm_session.picking_selections.picking_selections_list != [None,None,None,None]:
                        info = None
                        menu_type = "pick_menu"
                    else:
                        # Here the obj_menu is activated based on the atom that
                        # was identified by the picking function.
                        # The picking function detects the selected pixel and
                        # associates it with the respective object.
                        # There is no selection (blue dots) but an atom was
                        # identified in the click with the right button
                        if self.atom_picked is not None:
                            # Getting the info about the atom that was identified in the click
                            label = "{} / {} / {}({}) / {}({} / {})".format(self.atom_picked.vm_object.name,
                                                                        self.atom_picked.chain,
                                                                        self.atom_picked.resn,
                                                                        self.atom_picked.resi,
                                                                        self.atom_picked.name,
                                                                        self.atom_picked.index,
                                                                        self.atom_picked.symbol)
                            self.atom_picked = None
                            menu_type = "obj_menu"
                            info = label
                        else:
                            # When no atom is identified in the click (user clicked
                            # on a point in the background)
                            menu_type = "bg_menu"
                            info = None
                else:
                    # When a selection (viewing selection) is active, the
                    # selection menu is passed as an option.
                    info = self.vm_session.selections[self.vm_session.current_selection].get_selection_info()
                    menu_type = "sele_menu"
                # The right button (button = 3) always opens one of the available menus.
                self.parent_widget.show_gl_menu(menu_type=menu_type, info=info)
    
    def mouse_motion(self, mouse_x, mouse_y):
        """ Function doc
        """
        x = np.float32(mouse_x)
        y = np.float32(mouse_y)
        dx = x - self.mouse_x
        dy = y - self.mouse_y
        if (dx == 0) and (dy == 0):
            return
        self.mouse_x, self.mouse_y = x, y
        changed = False
        if self.mouse_rotate:
            changed = self._rotate_view(dx, dy, x, y)
        elif self.mouse_pan:
            changed = self._pan_view(x, y)
        elif self.mouse_zoom:
            changed = self._zoom_view(dy)
        if changed:
            self.dragging = True
            self.queue_draw()
    
    def mouse_scroll(self, direction):
        """ Function doc
        """
        up = np.int32(direction) == 1
        down = np.int32(direction) == -1
        if self.ctrl:
            if self.editing_mols:
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    if vm_object.editing:
                        if up:
                            vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, -self.scroll]))
                        if down:
                            vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, self.scroll]))
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        if vm_object.editing:
                            if up:
                                vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, -self.scroll]))
                            if down:
                                vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, self.scroll]))
            else:
                if up:
                    self.model_mat = mop.my_glTranslatef(self.model_mat, np.array([0.0, 0.0, -self.scroll]))
                if down:
                    self.model_mat = mop.my_glTranslatef(self.model_mat, np.array([0.0, 0.0, self.scroll]))
                
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    if up:
                        vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, -self.scroll]))
                    if down:
                        vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, self.scroll]))
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        if up:
                            vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, -self.scroll]))
                        if down:
                            vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, self.scroll]))
        else:
            pos_z = self.glcamera.get_position()[2]
            if up:
                self.glcamera.z_near -= self.scroll
                self.glcamera.z_far += self.scroll
            if down:
                if (self.glcamera.z_far-self.scroll) >= (self.glcamera.min_zfar):
                    if (self.glcamera.z_far-self.scroll) > (self.glcamera.z_near + self.scroll + 0.005):
                        self.glcamera.z_near += self.scroll
                        self.glcamera.z_far -= self.scroll
            
            if (self.glcamera.z_near >= self.glcamera.min_znear):
                self.glcamera.set_projection_matrix(mop.my_glPerspectivef(self.glcamera.field_of_view,
                        self.glcamera.viewport_aspect_ratio, self.glcamera.z_near, self.glcamera.z_far))
            else:
                if self.glcamera.z_far < (self.glcamera.min_zfar + self.glcamera.min_znear):
                    self.glcamera.z_near -= self.scroll
                    self.glcamera.z_far = self.glcamera.min_clip + self.glcamera.min_znear
                self.glcamera.set_projection_matrix(mop.my_glPerspectivef(self.glcamera.field_of_view,
                                                    self.glcamera.viewport_aspect_ratio,
                                                    self.glcamera.min_znear, self.glcamera.z_far))
            self.glcamera.update_fog()
        self.queue_draw()
    
    def _rotate_view(self, dx, dy, x, y):
        """ Function doc """
        angle = np.sqrt(dx**2 + dy**2) / (self.width + 1) * 180.0
        if self.shift:
            if self.show_dynamic_line:
                self.dynamic_line.end = self.get_viewport_pos(self.mouse_x, self.mouse_y)
                self.dynamic_line.update_points()
            else:
                self.selection_box.end = self.get_viewport_pos(self.mouse_x, self.mouse_y)
                self.selection_box.update_points()
        
        else:
            if self.ctrl:
                if abs(dx) >= abs(dy):
                    if (y - self.height / 2.0) < 0:
                        rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([0.0, 0.0, dx]))
                    else:
                        rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([0.0, 0.0, -dx]))
                else:
                    if (x - self.width / 2.0) < 0:
                        rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([0.0, 0.0, -dy]))
                    else:
                        rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([0.0, 0.0, dy]))
            else:
                rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([-dy, -dx, 0.0]))
            
            if self.editing_mols:
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    if vm_object.editing:
                        vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, rot_mat)
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        if vm_object.editing:
                            vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, rot_mat)
            else:
                self.model_mat = mop.my_glMultiplyMatricesf(self.model_mat, rot_mat)
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, rot_mat)
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, rot_mat)
            
            # Axis operations, this code only affects the gizmo axis
            if not self.editing_mols:
                self.axis.model_mat = mop.my_glTranslatef(self.axis.model_mat, -self.axis.zrp)
                if self.ctrl:
                    if abs(dx) >= abs(dy):
                        if (y - self.height / 2.0) < 0:
                            self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([0.0, 0.0, dx]))
                        else:
                            self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([0.0, 0.0, -dx]))
                    else:
                        if (x - self.width / 2.0) < 0:
                            self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([0.0, 0.0, -dy]))
                        else:
                            self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([0.0, 0.0, dy]))
                else:
                    self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([dy, dx, 0.0]))
                self.axis.model_mat = mop.my_glTranslatef(self.axis.model_mat, self.axis.zrp)
            # Axis operations, this code only affects the gizmo axis
        return True
    
    def _pan_view(self, x, y):
        """ Function doc """
        px, py, pz = self._mouse_pos(x, y)
        pan_mat = mop.my_glTranslatef(np.identity(4, dtype=np.float32), np.array(
                                    [(px - self.drag_pos_x) * self.glcamera.z_far / 10.0,
                                     (py - self.drag_pos_y) * self.glcamera.z_far / 10.0,
                                     (pz - self.drag_pos_z) * self.glcamera.z_far / 10.0]))
        if self.editing_mols:
            for index, vm_object in self.vm_session.vm_objects_dic.items():
                if vm_object.editing:
                    vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, pan_mat)
            
            for key in self.vm_session.vm_geometric_object_dic.keys():
                vm_object = self.vm_session.vm_geometric_object_dic[key]
                if vm_object:
                    if vm_object.editing:
                        vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, pan_mat)
        
        else:
            self.model_mat = mop.my_glMultiplyMatricesf(self.model_mat, pan_mat)
            for index, vm_object in self.vm_session.vm_objects_dic.items():
                vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, pan_mat)
            
            for key in self.vm_session.vm_geometric_object_dic.keys():
                vm_object = self.vm_session.vm_geometric_object_dic[key]
                if vm_object:
                    vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, pan_mat)
            self.zero_reference_point = mop.get_xyz_coords(self.model_mat)
        
        self.drag_pos_x = px
        self.drag_pos_y = py
        self.drag_pos_z = pz
        return True
    
    def _zoom_view(self, dy):
        """ Function doc """
        delta = (((self.glcamera.z_far - self.glcamera.z_near) / 2.0) + self.glcamera.z_near) / 200.0
        move_z = dy * delta
        moved_mat = mop.my_glTranslatef(self.glcamera.view_matrix, np.array([0.0, 0.0, move_z]))
        moved_pos = mop.get_xyz_coords(moved_mat)
        if moved_pos[2] > 0.101:
            self.glcamera.set_view_matrix(moved_mat)
            self.glcamera.z_near -= move_z
            self.glcamera.z_far -= move_z
            if self.glcamera.z_near >= self.glcamera.min_znear:
                self.glcamera.set_projection_matrix(mop.my_glPerspectivef(self.glcamera.field_of_view, 
                                                    self.glcamera.viewport_aspect_ratio,
                                                    self.glcamera.z_near, self.glcamera.z_far))
            else:
                if self.glcamera.z_far < (self.glcamera.min_zfar+self.glcamera.min_znear):
                    self.glcamera.z_near += move_z
                    self.glcamera.z_far = self.glcamera.min_zfar+self.glcamera.min_znear
                self.glcamera.set_projection_matrix(mop.my_glPerspectivef(self.glcamera.field_of_view, 
                                                    self.glcamera.viewport_aspect_ratio,
                                                    self.glcamera.min_znear, self.glcamera.z_far))
            self.glcamera.update_fog()
            self.dist_cam_zrp += -move_z
            return True
        return False
    
    def render(self):
        """ This is the function that will be called everytime the window
            needs to be re-drawed.
        """
        if self.shader_flag:
            self.create_gl_programs()
            self.selection_box.initialize_gl()
            self.axis.initialize_gl()
            self.shader_flag = False
        if self.selection_box_picking:
            self._selection_box_pick()
        if self.picking:
            self._pick()
        
        GL.glClearColor(self.bckgrnd_color[0], self.bckgrnd_color[1],
                        self.bckgrnd_color[2], self.bckgrnd_color[3])
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        
        for index, vm_object in self.vm_session.vm_objects_dic.items():
            if vm_object.active:
                if len(vm_object.frames) > 0:
                    for rep_name in vm_object.representations:
                        if vm_object.representations[rep_name] is not None:
                            # Only shows the representation if
                            # representations[rep_name].active = True
                            if vm_object.representations[rep_name].active:
                                print(rep_name, "<- Debugging")
                                vm_object.representations[rep_name].draw_representation()
        # Check if the picking function is active.
        # Viewing and picking selections cannot be displayed at the same time.
        if self.vm_session._picking_selection_mode:
            self._draw_picking_label()
            for rep_name in self.vm_session.vm_geometric_object_dic.keys():
                if self.vm_session.vm_geometric_object_dic[rep_name]:
                    if self.vm_session.vm_geometric_object_dic[rep_name].representations["dotted_lines"].active:
                        self.vm_session.vm_geometric_object_dic[rep_name].representations["dotted_lines"].draw_representation()
        else:
            for vm_object in self.vm_session.selections[self.vm_session.current_selection].selected_objects:
                # Here are represented the blue dots referring to the atom's selections
                if vm_object.selection_dots_vao is None:
                    shapes._make_gl_selection_dots(self.shader_programs["picking_dots"],
                                                   vismol_object=vm_object)
                
                # Extracting the indexes for each vismol_object that was selected
                indexes = self.vm_session.selections[self.vm_session.current_selection].selected_objects[vm_object]
                size = self.vm_config.gl_parameters["dot_sel_size"]
               
                GL.glPointSize(size * self.height / (abs(self.dist_cam_zrp)) / 2)
                GL.glUseProgram(self.shader_programs["picking_dots"])
                GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
                self.load_matrices(self.shader_programs["picking_dots"], vm_object.model_mat)
                GL.glBindVertexArray(vm_object.selection_dots_vao)
                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vm_object.selection_dot_buffers[0])
                GL.glBufferData(GL.GL_ARRAY_BUFFER, indexes.itemsize * len(indexes),
                                indexes, GL.GL_STATIC_DRAW)
                frame = self._safe_frame_exchange(vm_object)
                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vm_object.selection_dot_buffers[1])
                GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.itemsize * len(frame),
                                frame, GL.GL_STATIC_DRAW)
                GL.glDrawElements(GL.GL_POINTS, np.int32(len(indexes)), GL.GL_UNSIGNED_INT, None)
                GL.glBindVertexArray(0)
                GL.glDisable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
                GL.glPointSize(1)
                GL.glUseProgram(0)
                GL.glDisable(GL.GL_DEPTH_TEST)
        
        if self.show_dynamic_line and self.shift:
            if self.dynamic_line.vao is None:
                self.dynamic_line._make_gl_selection_box()
            else:
                self.dynamic_line._draw_selection_box()
        
        if self.show_selection_box and self.shift:
            if self.selection_box.vao is None:
                self.selection_box._make_gl_selection_box()
            else:
                self.selection_box._draw()
        
        if self.show_axis:
            self.axis._draw(True)
            self.axis._draw(False)
        return True
    
    def create_gl_programs(self):
        """ Function doc
        """
        print("OpenGL version: ", GL.glGetString(GL.GL_VERSION))
        print("OpenGL major version: ", GL.glGetDoublev(GL.GL_MAJOR_VERSION))
        print("OpenGL minor version: ", GL.glGetDoublev(GL.GL_MINOR_VERSION))
        #----------------------------------------------------------------------#
        for rep in self.representations_available:
            func = getattr(self, "_compile_shader_" + rep)
            try:
                func()
            except AttributeError as ae:
                print("Representation of type '{}' not implemented".format(rep))
                print(ae)
    
    def load_shaders(self, vertex, fragment, geometry=None):
        """ Here the shaders are loaded and compiled to an OpenGL program. By default
            the constructor shaders will be used, if you want to change the shaders
            use this function. The flag is used to create only one OpenGL program.
            
            Keyword arguments:
            vertex -- The vertex shader to be used
            fragment -- The fragment shader to be used
        """
        my_vertex_shader = self.create_shader(vertex, GL.GL_VERTEX_SHADER)
        my_fragment_shader = self.create_shader(fragment, GL.GL_FRAGMENT_SHADER)
        if geometry is not None:
            my_geometry_shader = self.create_shader(geometry, GL.GL_GEOMETRY_SHADER)
        program = GL.glCreateProgram()
        GL.glAttachShader(program, my_vertex_shader)
        GL.glAttachShader(program, my_fragment_shader)
        if geometry is not None:
            GL.glAttachShader(program, my_geometry_shader)
        GL.glLinkProgram(program)
        return program
    
    def create_shader(self, shader_prog, shader_type):
        """ Creates, links to a source, compiles and returns a shader.
            
            Keyword arguments:
            shader -- The shader text to use
            shader_type -- The OpenGL enum type of shader, it can be:
                           GL.GL_VERTEX_SHADER, GL.GL_GEOMETRY_SHADER or GL.GL_FRAGMENT_SHADER
            
            Returns:
            A shader object identifier or pops out an error
        """
        shader = GL.glCreateShader(shader_type)
        GL.glShaderSource(shader, shader_prog)
        GL.glCompileShader(shader)
        if GL.glGetShaderiv(shader, GL.GL_COMPILE_STATUS) != GL.GL_TRUE:
            print("Error compiling the shader: ", shader_type)
            raise RuntimeError(GL.glGetShaderInfoLog(shader))
        return shader
    
    def _selection_box_pick(self):
        """ Selects a set of atoms from pixels obtained by the selection rectangle.  
            This function (method) is called in the render method, when the 
            "self.selection_box_picking" attribute is active. 
        
        """
        
        # glReadPixels and glReadnPixels return pixel data from the frame buffer, 
        # starting with the pixel whose lower left corner is at location (x, y), 
        # into client memory starting at location data.
        #
        # In GTK, x=0 and y=0 set to upper left corner (unlike openGL input data, 
        # the following lines do the coordinate conversion) 
        
        selection_box_x2 = self.mouse_x
        selection_box_y2 = self.height - self.mouse_y
        selection_box_width  = selection_box_x2 - self.selection_box_x
        selection_box_height = selection_box_y2 - self.selection_box_y
        
        #Looking for the lower left corner of the checkbox
        if selection_box_width > 0 and selection_box_height > 0:
            pos_x = self.selection_box_x
            pos_y = self.selection_box_y
            width = selection_box_width
            height = selection_box_height
        
        elif selection_box_width < 0 and selection_box_height > 0:
            pos_x = selection_box_x2
            pos_y = self.selection_box_y
            width = -selection_box_width
            height = selection_box_height
        
        elif selection_box_width < 0 and selection_box_height < 0:
            pos_x = selection_box_x2
            pos_y = selection_box_y2
            width =  -selection_box_width
            height = -selection_box_height
        else:
            pos_x = self.selection_box_x
            pos_y = selection_box_y2
            width = selection_box_width
            height = -selection_box_height
        
        # taking the module from the width and height values 
        if pos_x < 0:
            pos_x = 0.0
        if pos_y < 0:
            pos_y = 0.0
        
        GL.glClearColor(1, 1, 1, 1)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        for index, vm_object in self.vm_session.vm_objects_dic.items():
            if vm_object.active:
                #vismol_object has few different types of representations
                for rep_name in vm_object.representations:
                    # checking all the representations in vismol_object.representations dictionary
                    if vm_object.representations[rep_name] is not None:
                        #  vismol_object.representations[rep_name] may be active or not  True/False
                        if vm_object.representations[rep_name].active:
                            vm_object.representations[rep_name].draw_background_sel_representation()
        
        GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
        data = GL.glReadPixels(pos_x, pos_y, width, height, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
        data = list(data)
        picked_set = set()
        for i in range(0, len(data), 4):
            #converting RGB values to atoms address (unique id)
            pickedID = data[i] + data[i+1] * 256 + data[i+2] * 256 * 256;
            picked_set.add(pickedID)
        for pickedID in picked_set:
            if pickedID == 16777215:
                pass
            else:
                self.atom_picked = self.vm_session.atom_dic_id[pickedID]
                # The disable variable does not allow, if the selected 
                # atom is already in the selected list, to be removed.
                # The disable variable is "False" for when we use 
                # selection by area (selection box)
                self.vm_session._selection_function(selected=self.atom_picked, disable=False)
        self.selection_box_picking = False
    
    def _pick(self):
        """ Function doc """
        GL.glClearColor(1, 1, 1, 1)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        for index, vm_object in self.vm_session.vm_objects_dic.items():
            if vm_object.active:
                #vismol_object has few different types of representations
                for rep_name in vm_object.representations:
                    # checking all the representations in vismol_object.representations dictionary
                    if vm_object.representations[rep_name] is None:
                        pass
                    else:
                        #  vismol_object.representations[rep_name] may be active or not  True/False
                        if vm_object.representations[rep_name].active:
                            vm_object.representations[rep_name].draw_background_sel_representation()
        
        GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
        pos = [self.picking_x, self.height - self.picking_y]
        data = GL.glReadPixels(pos[0], pos[1], 1, 1, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
        
        #converting RGB values to atoms address (unique id)
        pickedID = data[0] + data[1] * 256 + data[2] * 256 * 256;
        if pickedID == 16777215:
            self.atom_picked = None
            if self.button == 1:
                self.vm_session._selection_function(self.atom_picked)
                self.button = None
        else:
            try:
                # Using antialias, in some rare cases, the pick function is not 
                # identifying the right color of the selected atom. This event is 
                # rare, but can impair viewing if it is not properly ignored
                self.atom_picked = self.vm_session.atom_dic_id[pickedID]
                if self.button == 1:
                    self.vm_session._selection_function(self.atom_picked)
                    self.button = None
            except KeyError as ke:
                print("pickedID", pickedID, "not found")
                print(ke)
                self.button = None
        self.picking = False
        return True
    
    def load_fog(self, program):
        """ Load the fog parameters in the specified program
            
            fog_start -- The coordinates where the fog will begin (always
                         positive)
            fog_end -- The coordinates where the fog will begin (always positive
                       and greater than fog_start)
            fog_color -- The color for the fog (same as background)
        """
        fog_s = GL.glGetUniformLocation(program, "fog_start")
        GL.glUniform1fv(fog_s, 1, self.glcamera.fog_start)
        fog_e = GL.glGetUniformLocation(program, "fog_end")
        GL.glUniform1fv(fog_e, 1, self.glcamera.fog_end)
        fog_c = GL.glGetUniformLocation(program, "fog_color")
        GL.glUniform4fv(fog_c, 1, self.bckgrnd_color)
    
    def load_matrices(self, program=None, model_mat=None):
        """ Load the matrices to OpenGL.
            
            model_mat -- transformation matrix for the objects rendered
            view_mat -- transformation matrix for the camera used
            proj_mat -- matrix for the space to be visualized in the scene
        """
        model = GL.glGetUniformLocation(program, "model_mat")
        GL.glUniformMatrix4fv(model, 1, GL.GL_FALSE, model_mat)
        view = GL.glGetUniformLocation(program, "view_mat")
        GL.glUniformMatrix4fv(view, 1, GL.GL_FALSE, self.glcamera.view_matrix)
        proj = GL.glGetUniformLocation(program, "proj_mat")
        GL.glUniformMatrix4fv(proj, 1, GL.GL_FALSE, self.glcamera.projection_matrix)
    
    def load_dot_params(self, program):
        """ Function doc
        """
        # Extern line
        linewidth = np.float32(80 / abs(self.dist_cam_zrp))
        if linewidth > 3.73:
            linewidth = 3.73
        # Intern line
        antialias = np.float32(80 / abs(self.dist_cam_zrp))
        if antialias > 3.73:
            antialias = 3.73
        # Dot size factor
        dot_factor = np.float32(500 / abs(self.dist_cam_zrp))
        if dot_factor > 150.0:
            dot_factor = 150.0
        uni_vext_linewidth = GL.glGetUniformLocation(program, "vert_ext_linewidth")
        GL.glUniform1fv(uni_vext_linewidth, 1, linewidth)
        uni_vint_antialias = GL.glGetUniformLocation(program, "vert_int_antialias")
        GL.glUniform1fv(uni_vint_antialias, 1, antialias)
        uni_dot_size = GL.glGetUniformLocation(program, "vert_dot_factor")
        GL.glUniform1fv(uni_dot_size, 1, dot_factor)
    
    def load_lights(self, program):
        """ Function doc
        """
        light_pos = GL.glGetUniformLocation(program, "my_light.position")
        GL.glUniform3fv(light_pos, 1, self.light_position)
        amb_coef = GL.glGetUniformLocation(program, "my_light.ambient_coef")
        GL.glUniform1fv(amb_coef, 1, self.light_ambient_coef)
        shiny = GL.glGetUniformLocation(program, "my_light.shininess")
        GL.glUniform1fv(shiny, 1, self.light_shininess)
        intensity = GL.glGetUniformLocation(program, "my_light.intensity")
        GL.glUniform3fv(intensity, 1, self.light_intensity)
    
    def load_antialias_params(self, program):
        """ Function doc """
        a_length = GL.glGetUniformLocation(program, "antialias_length")
        GL.glUniform1fv(a_length, 1, 0.05)
        bck_col = GL.glGetUniformLocation(program, "alias_color")
        GL.glUniform3fv(bck_col, 1, self.bckgrnd_color[:3])
    
    def _draw_picking_label(self):
        """ This function draws the labels of the atoms selected by the
            function picking #1 #2 #3 #4
        """
        if self.vm_font.vao is None:
            self.vm_font.make_freetype_font()
            self.vm_font.make_freetype_texture(self.shader_programs["freetype"])
        
        number = 1
        self.chars     = 0
        self.xyz_pos   = []
        self.uv_coords = []
        
        for atom in self.vm_session.picking_selections.picking_selections_list:
            if atom:
                text = "#" + str(number)
                frame = self._get_vismol_object_frame(atom.vm_object)
                x, y, z = atom.coords(frame)
                point = np.array([x, y, z, 1], dtype=np.float32)
                point = np.dot(point, self.model_mat)
                GL.glBindTexture(GL.GL_TEXTURE_2D, self.vm_font.texture_id)
                for i, c in enumerate(text):
                    self.chars += 1
                    c_id = ord(c)
                    x = c_id % 16
                    y = c_id // 16 - 2
                    self.xyz_pos.append(point[0] + i * self.vm_font.char_width)
                    self.xyz_pos.append(point[1])
                    self.xyz_pos.append(point[2])
                    self.uv_coords.append(x * self.vm_font.text_u)
                    self.uv_coords.append(y * self.vm_font.text_v)
                    self.uv_coords.append((x + 1) * self.vm_font.text_u)
                    self.uv_coords.append((y + 1) * self.vm_font.text_v)
            number += 1
        self.xyz_pos = np.array(self.xyz_pos, dtype=np.float32)
        self.uv_coords = np.array(self.uv_coords, dtype=np.float32)
        
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vm_font.vbos[0])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.xyz_pos.itemsize * len(self.xyz_pos),
                        self.xyz_pos, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vm_font.vbos[1])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.uv_coords.itemsize * len(self.uv_coords),
                        self.uv_coords, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glUseProgram(self.shader_programs["freetype"])
        
        self.vm_font.load_matrices(self.shader_programs["freetype"],
                                   self.glcamera.view_matrix,
                                   self.glcamera.projection_matrix)
        self.vm_font.load_font_params(self.shader_programs["freetype"])
        
        GL.glBindVertexArray(self.vm_font.vao)
        GL.glDrawArrays(GL.GL_POINTS, 0, self.chars)
        GL.glDisable(GL.GL_BLEND)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
    
    def _compile_shader_picking_dots(self):
        """ Function doc """
        self.shader_programs["picking_dots"] = self.load_shaders(shaders_pick.vertex_shader_picking_dots,
                                                                 shaders_pick.fragment_shader_picking_dots)
    
    def _compile_shader_dots(self):
        """ Function doc """
        dot_type = self.vm_config.gl_parameters["dot_type"]
        self.shader_programs["dots"] = self.load_shaders(shaders_dots.shader_type[dot_type]["vertex_shader"],
                                                         shaders_dots.shader_type[dot_type]["fragment_shader"])
        self.shader_programs["dots_sel"] = self.load_shaders(shaders_dots.shader_type[dot_type]["sel_vertex_shader"],
                                                             shaders_dots.shader_type[dot_type]["sel_fragment_shader"])
    
    def _compile_shader_lines(self):
        """ Function doc """
        line_type = self.vm_config.gl_parameters["line_type"]
        self.shader_programs["lines"] = self.load_shaders(shaders_lines.shader_type[line_type]["vertex_shader"],
                                                          shaders_lines.shader_type[line_type]["fragment_shader"],
                                                          shaders_lines.shader_type[line_type]["geometry_shader"])
        self.shader_programs["lines_sel"] = self.load_shaders(shaders_lines.shader_type[line_type]["sel_vertex_shader"],
                                                              shaders_lines.shader_type[line_type]["sel_fragment_shader"],
                                                              shaders_lines.shader_type[line_type]["sel_geometry_shader"])
    
    def _compile_shader_nonbonded(self):
        """ Function doc """
        self.shader_programs["nonbonded"] = self.load_shaders(shaders_nonbonded.vertex_shader_non_bonded,
                                                              shaders_nonbonded.fragment_shader_non_bonded,
                                                              shaders_nonbonded.geometry_shader_non_bonded)
        self.shader_programs["nonbonded_sel"] = self.load_shaders(shaders_nonbonded.sel_vertex_shader_non_bonded,
                                                                  shaders_nonbonded.sel_fragment_shader_non_bonded,
                                                                  shaders_nonbonded.sel_geometry_shader_non_bonded)
    
    def _compile_shader_dotted_lines(self):
        """ Function doc """
        line_type = 3
        self.shader_programs["dotted_lines"] = self.load_shaders(shaders_lines.shader_type[line_type]["vertex_shader"],
                                                                 shaders_lines.shader_type[line_type]["fragment_shader"],
                                                                 shaders_lines.shader_type[line_type]["geometry_shader"])
        self.shader_programs["dotted_lines_sel"] = self.load_shaders( shaders_lines.shader_type[line_type]["sel_vertex_shader"],
                                                                      shaders_lines.shader_type[line_type]["sel_fragment_shader"],
                                                                      shaders_lines.shader_type[line_type]["sel_geometry_shader"])
    
    def _compile_shader_ribbon(self):
        """ Function doc """
        line_type = self.vm_config.gl_parameters["ribbon_type"]
        self.shader_programs["ribbon"] = self.load_shaders(shaders_lines.shader_type[line_type]["vertex_shader"],
                                                            shaders_lines.shader_type[line_type]["fragment_shader"],
                                                            shaders_lines.shader_type[line_type]["geometry_shader"])
        self.shader_programs["ribbon_sel"] = self.load_shaders(shaders_lines.shader_type[line_type]["sel_vertex_shader"],
                                                                shaders_lines.shader_type[line_type]["sel_fragment_shader"],
                                                                shaders_lines.shader_type[line_type]["sel_geometry_shader"])
    
    def _compile_shader_sticks(self):
        """ Function doc """
        self.shader_programs["sticks"] = self.load_shaders(shaders_sticks.vertex_shader_sticks,
                                                           shaders_sticks.fragment_shader_sticks,
                                                           shaders_sticks.geometry_shader_sticks)
        self.shader_programs["sticks_sel"] = self.load_shaders(shaders_sticks.sel_vertex_shader_sticks,
                                                               shaders_sticks.sel_fragment_shader_sticks,
                                                               shaders_sticks.sel_geometry_shader_sticks)
    
    def _compile_shader_spheres(self):
        """ Function doc """
        self.shader_programs["spheres"] = self.load_shaders(shaders_spheres.vertex_shader_spheres,
                                                            shaders_spheres.fragment_shader_spheres)
        self.shader_programs["spheres_sel"] = self.load_shaders(shaders_spheres.vertex_shader_spheres,
                                                                shaders_spheres.fragment_shader_spheres)
    
    def _compile_shader_impostor(self):
        """ Function doc """
        im_type = 2
        self.shader_programs["impostor"] = self.load_shaders(shaders_impostor.shader_type[im_type]["vertex_shader"],
                                                             shaders_impostor.shader_type[im_type]["fragment_shader"],
                                                             shaders_impostor.shader_type[im_type]["geometry_shader"])
        self.shader_programs["impostor_sel"] = self.load_shaders(shaders_impostor.shader_type[0]["sel_vertex_shader"],
                                                               shaders_impostor.shader_type[0]["sel_fragment_shader"])
    
    def _compile_shader_surface(self):
        """ Function doc """
        self.shader_programs["surface"] = self.load_shaders(shaders_surface.vertex_shader_surface,
                                                            shaders_surface.fragment_shader_surface,
                                                            shaders_surface.geometry_shader_surface)
        self.shader_programs["surface_sel"] = self.load_shaders(shaders_spheres.vertex_shader_spheres,
                                                                shaders_spheres.fragment_shader_spheres)
        
    def _compile_shader_cartoon(self):
        """ Function doc """
        self.shader_programs["cartoon"] = self.load_shaders(shaders_cartoon.v_shader_triangles,
                                                            shaders_cartoon.f_shader_triangles)
    
    def _compile_shader_freetype(self):
        """ Function doc """
        self.shader_programs["freetype"] = self.load_shaders(shaders_vm_freetype.vertex_shader_freetype,
                                                             shaders_vm_freetype.fragment_shader_freetype,
                                                             shaders_vm_freetype.geometry_shader_freetype)
    
    #----------------------------NOT IMPLEMENTED YET---------------------------#
    def _dynamic_bonds_shaders(self):
        """ Function doc """
        self.shader_programs["dynamic"] = self.load_shaders(shaders_sticks.vertex_shader_sticks,
                                                            shaders_sticks.fragment_shader_sticks,
                                                            shaders_sticks.geometry_shader_sticks)
        self.shader_programs["dynamic_sel"] = self.load_shaders(shaders_sticks.sel_vertex_shader_sticks,
                                                               shaders_sticks.sel_fragment_shader_sticks,
                                                               shaders_sticks.sel_geometry_shader_sticks)
    
    def _wires_dot_shaders(self):
        """ Function doc """
        self.shader_programs["wires"] = self.load_shaders(shaders_wires.vertex_shader_wires,
                                                          shaders_wires.fragment_shader_wires,
                                                          shaders_wires.geometry_shader_wires)
        self.shader_programs["wires_sel"] = self.load_shaders(shaders_spheres.vertex_shader_spheres,
                                                              shaders_spheres.fragment_shader_spheres)
    #----------------------------NOT IMPLEMENTED YET---------------------------#
    
    def _safe_frame_exchange(self, vismol_object):
        """ Function doc 
        This function checks if the number of the called frame will not exceed 
        the limit of frames that each object has. Allowing two objects with 
        different trajectory sizes to be manipulated at the same time within the 
        glArea
        """
        if self.frame < 0:
            self.frame = 0
        if self.frame >= (len(vismol_object.frames) - 1):
            frame = vismol_object.frames[len(vismol_object.frames) - 1]
        else:
            frame = vismol_object.frames[self.frame]
        return frame
    
    def _get_vismol_object_frame(self, vismol_object):
        """ Function doc """
        if self.frame < 0:
            self.frame = 0
        if self.frame >= (len(vismol_object.frames) - 1):
            frame = len(vismol_object.frames) - 1
        else:
            frame = self.frame
        return frame
    
    def get_viewport_pos(self, x, y):
        """ Function doc """
        px = (2.0 * x - self.width) / self.width
        py = (2.0 * y - self.height) / self.height
        return np.array([px, -py], dtype=np.float32)
    
    def _mouse_pos(self, x, y):
        """
        Use the ortho projection and viewport information
        to map from mouse co-ordinates back into world
        co-ordinates
        """
        px = x / self.width
        py = y / self.height
        px = self.left + px * (self.right - self.left)
        py = self.top + py * (self.bottom - self.top)
        pz = self.glcamera.z_near
        return px, py, pz
    
    def center_on_atom(self, atom):
        """ Function doc
        """
        frame_index = self._get_vismol_object_frame(atom.vm_object)
        self.center_on_coordinates(atom.vm_object, atom.coords(frame_index))
        return True
    
    def center_on_coordinates(self, vismol_object, target):
        """ Takes the coordinates of an atom in absolute coordinates and first
            transforms them in 4D world coordinates, then takes the unit vector
            of that atom position to generate the loop animation. To generate
            the animation, first obtains the distance from the zero reference
            point (always 0,0,0) to the atom, then divides this distance in a
            defined number of cycles, this result will be the step for
            translation. For the translation, the world will move a number of
            steps defined, and every new point will be finded by multiplying the
            unit vector by the step. As a final step, to avoid biases, the world
            will be translated to the atom position in world coordinates.
            The effects will be applied on the model matrices of every VisMol
            object and the model matrix of the window.
        """
        if (self.zero_reference_point[0] != target[0]) or \
           (self.zero_reference_point[1] != target[1]) or \
           (self.zero_reference_point[2] != target[2]):
            self.zero_reference_point[:] = target
            pos = np.array([target[0],target[1],target[2],1], dtype=np.float32)
            model_pos = vismol_object.model_mat.T.dot(pos)[:3]
            self.model_mat = mop.my_glTranslatef(self.model_mat, -model_pos)
            unit_vec = model_pos / np.linalg.norm(model_pos)
            step = np.linalg.norm(model_pos)/15.0
            for i in range(15):
                to_move = unit_vec * step
                
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -to_move)
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -to_move)
                
                # WARNING: Method only works with GTK!!!
                self.parent_widget.get_window().invalidate_rect(None, False)
                self.parent_widget.get_window().process_updates(False)
                # WARNING: Method only works with GTK!!!
                time.sleep(self.vm_config.gl_parameters["center_on_coord_sleep_time"])
            
            for index, vm_object in self.vm_session.vm_objects_dic.items():
                model_pos = vm_object.model_mat.T.dot(pos)[:3]
                vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -model_pos)
            
            for key in self.vm_session.vm_geometric_object_dic.keys():
                vm_object = self.vm_session.vm_geometric_object_dic[key]
                if vm_object:
                    model_pos = vm_object.model_mat.T.dot(pos)[:3]
                    vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -model_pos)
            
            self.queue_draw()
        return True
    
    def queue_draw(self):
        """ Function doc """
        self.parent_widget.queue_draw()
    
    def _print_matrices(self):
        """ Function doc
        """
        print(self.model_mat,"<== widget model_mat")
        for index, vm_object in self.vm_session.vm_objects_dic.items():
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print(vm_object.model_mat,"<== vismol_object model_mat")
    
