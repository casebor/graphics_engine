#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from OpenGL import GL
from logging import getLogger
from vismol.libgl.vismol_gl_core import VismolGLCore
import vismol.utils.matrix_operations as mop
import vismol.libgl.shaders.sticks as shaders_sticks
import vismol.libgl.shaders.spheres as shaders_spheres


logger = getLogger(__name__)


class VismolGLBuilder(VismolGLCore):
    
    def __init__(self, widget, vismol_config, width=640.0, height=420.0):
        super(VismolGLBuilder, self).__init__(widget, vismol_config, width, height)
        self.vm_session = None
    
    def initialize_builder(self):
        self.initialize_parent()
    
    def _compile_shader_sticks(self):
        """ Function doc """
        sticks_type = self.vm_config.gl_parameters["sticks_type"]
        self.shader_programs["sticks"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["vertex_shader"],
                                                   shaders_sticks.shader_type[sticks_type]["fragment_shader"],
                                                   shaders_sticks.shader_type[sticks_type]["geometry_shader"])
        self.shader_programs["sticks_sel"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["sel_vertex_shader"],
                                                       shaders_sticks.shader_type[sticks_type]["sel_fragment_shader"],
                                                       shaders_sticks.shader_type[sticks_type]["sel_geometry_shader"])
    
    def _compile_shader_spheres(self):
        """ Function doc """
        self.shader_programs["spheres"] = self.load_shaders(shaders_spheres.vertex_shader_spheres,
                                                    shaders_spheres.fragment_shader_spheres)
        self.shader_programs["spheres_sel"] = self.load_shaders(shaders_spheres.sel_vertex_shader_spheres,
                                                        shaders_spheres.sel_fragment_shader_spheres)
    
    def mouse_pressed(self, button_number, mouse_x, mouse_y):
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
            # print("Left mouse pressed")
            pass
        if middle:
            # print("Middle mouse pressed")
            pass
        if right:
            # print("Right mouse pressed")
            pass
    
    def mouse_released(self, button_number, mouse_x, mouse_y):
        """ Function doc
        """
        left   = np.int32(button_number) == 1
        middle = np.int32(button_number) == 2
        right  = np.int32(button_number) == 3
        self.mouse_rotate = False
        self.mouse_zoom = False
        self.mouse_pan = False
        if self.dragging:
            pass
        else:
            if left:
                self.add_atom(mouse_x, mouse_y)
                self.mouse_button = 1
                #dragging is set to false here
                self.dragging = False
                self.parent_widget.queue_draw()
            if middle:
                pass
            if right:
                pass
        self.dragging = False
        self.parent_widget.queue_draw()
    
    def mouse_scroll(self, direction):
        """ Function doc
        """
        up = np.int32(direction) == 1
        down = np.int32(direction) == -1
        pos_z = self.glcamera.get_position()[2]
        if up:
            self.glcamera.z_near -= self.scroll
            self.glcamera.z_far += self.scroll
        if down:
            if (self.glcamera.z_far-self.scroll) >= self.glcamera.min_zfar:
                if (self.glcamera.z_far-self.scroll) > (self.glcamera.z_near + self.scroll + 0.005):
                    self.glcamera.z_near += self.scroll
                    self.glcamera.z_far -= self.scroll
        
        if self.glcamera.z_near >= self.glcamera.min_znear:
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
        self.parent_widget.queue_draw()
    
    def _rotate_view(self, dx, dy, x, y):
        """ Function doc """
        angle = np.sqrt(dx**2 + dy**2) / (self.width + 1) * 180.0
        rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([-dy, -dx, 0.0], dtype=np.float32))
        
        self.model_mat = mop.my_glMultiplyMatricesf(self.model_mat, rot_mat)
        for index, vm_object in self.vm_session.vm_objects_dic.items():
            vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, rot_mat)
        
        # Axis operations, this code only affects the gizmo axis
        self.axis.model_mat = mop.my_glTranslatef(self.axis.model_mat, -self.axis.zrp)
        self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([dy, dx, 0.0], dtype=np.float32))
        self.axis.model_mat = mop.my_glTranslatef(self.axis.model_mat, self.axis.zrp)
        # Axis operations, this code only affects the gizmo axis
        return True
    
    def _pan_view(self, x, y):
        px, py, pz = self._mouse_pos(x, y)
        pan_mat = mop.my_glTranslatef(np.identity(4, dtype=np.float32), np.array(
                                    [(px - self.drag_pos_x) * self.glcamera.z_far / 10.0,
                                     (py - self.drag_pos_y) * self.glcamera.z_far / 10.0,
                                     (pz - self.drag_pos_z) * self.glcamera.z_far / 10.0]))
        
        self.model_mat = mop.my_glMultiplyMatricesf(self.model_mat, pan_mat)
        for vm_object in self.vm_session.vm_objects_dic.values():
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
        if self.shader_flag:
            self.create_gl_programs()
            self.axis.initialize_gl()
            self.shader_flag = False
        GL.glClearColor(0.0, 0.0, 0.0, 1.0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        
        for vm_object in self.vm_session.vm_objects_dic.values():
            if vm_object.active:
                for representation in vm_object.representations.values():
                    if representation is not None:
                        if representation.active:
                            representation.draw_representation()
        
        if self.show_axis:
            self.axis._draw(True)
            self.axis._draw(False)
        return True
    
    def create_gl_programs(self):
        """ Function doc
        """
        for rep in self.vm_config.representations_available:
            func = getattr(self, "_compile_shader_" + rep)
            try:
                func()
            except AttributeError as ae:
                logger.error("Representation of type '%s' not implemented", rep)
                logger.error(ae)
    
    def add_atom(self, mouse_x, mouse_y):
        """ Function doc """
        cam_pos = self.glcamera.get_modelview_position(self.model_mat)
        proj = np.asmatrix(self.glcamera.projection_matrix)
        view = np.asmatrix(self.glcamera.view_matrix)
        model = np.asmatrix(self.model_mat)
        i_proj = proj.I
        i_view = view.I
        i_model = model.I
        i_mvp = i_proj * i_view * i_model
        viewport_pos = self.get_viewport_pos(mouse_x, mouse_y)
        viewport_pos.append(1)
        viewport_pos = np.asmatrix(viewport_pos)
        viewport_pos = (viewport_pos*i_mvp).A1
        viewport_pos /= viewport_pos[3]
        u_vec = mop.get_unit_vector(viewport_pos[:3] - cam_pos)
        v_vec = mop.get_unit_vector(-cam_pos)
        angle = mop.get_angle_rad(v_vec, u_vec)
        hypo = mop.get_euclidean(cam_pos, np.zeros(3, dtype=np.float32)) / np.cos(angle)
        xyz_pos = cam_pos + (u_vec * hypo)
        self.vm_session.add_new_atom(np.array(xyz_pos[:3], dtype=np.float32))
        self.queue_draw()
    
    def _safe_frame_coords(self, vm_object):
        """ Function doc
        """
        frame_coords = np.array(vm_object.molecule.frames[0])
        return frame_coords, 0
    
