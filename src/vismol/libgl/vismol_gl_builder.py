#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import time
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
    
    # def _compile_shader_sticks(self):
    #     """ Function doc """
    #     sticks_type = self.vm_config.gl_parameters["sticks_type"]
    #     self.shader_programs["sticks"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["vertex_shader"],
    #                                                shaders_sticks.shader_type[sticks_type]["fragment_shader"],
    #                                                shaders_sticks.shader_type[sticks_type]["geometry_shader"])
    #     self.shader_programs["sticks_sel"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["sel_vertex_shader"],
    #                                                    shaders_sticks.shader_type[sticks_type]["sel_fragment_shader"],
    #                                                    shaders_sticks.shader_type[sticks_type]["sel_geometry_shader"])
    
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
        pass
    
    def _zoom_view(self, dy):
        pass
    
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
        mod = self.get_viewport_pos(mouse_x, mouse_y)
        mod.append(1)
        mod = np.asmatrix(mod)
        mod = (mod*i_mvp).A1
        mod /= mod[3]
        u_vec = self.unit_vector(mod[:3] - cam_pos)
        v_vec = self.unit_vector(-cam_pos)
        angle = np.radians(self.get_angle(v_vec, u_vec))
        hypo = self.get_euclidean(cam_pos, [0,0,0]) / np.cos(angle)
        test = u_vec * hypo
        mod = cam_pos + test
        self.vm_session.add_new_atom(np.array(mod[:3], dtype=np.float32))
        self.queue_draw()
    
    def get_viewport_pos(self, x, y):
        """ Function doc """
        px = (2.0*x - self.width)/self.width
        py = (self.height - 2.0*y)/self.height
        return [px, py, self.glcamera.z_near]
    
    def unit_vector(self, vector):
        """ Returns the unit vector of the vector.
        """
        return vector / np.linalg.norm(vector)
    
    def get_angle(self, vecA, vecB):
        """ Return the angle in degrees of two vectors.
        """
        vecA_u = self.unit_vector(vecA)
        vecB_u = self.unit_vector(vecB)
        return np.degrees(np.arccos(np.clip(np.dot(vecA_u, vecB_u), -1.0, 1.0)))
    
    def get_euclidean(self, pa, pb):
        """ Returns the distance between two points in R3
        """
        import math
        if int(len(pa)) == 1:
            pa = [pa[0], 0.0, 0.0]
        if int(len(pa)) == 2:
            pa = [pa[0], pa[1], 0.0]
        if int(len(pb)) == 1:
            pb = [pb[0], 0.0, 0.0]
        if int(len(pb)) == 2:
            pb = [pb[0], pb[1], 0.0]
        return math.sqrt((pb[0]-pa[0])**2 + (pb[1]-pa[1])**2 + (pb[2]-pa[2])**2)
    
    def _safe_frame_coords(self, vm_object):
        """ Function doc
        """
        frame_coords = np.array(vm_object.molecule.frames[0])
        return frame_coords, 0
    
    # def center_on_coordinates(self, vismol_object, target):
    #     """ Takes the coordinates of an atom in absolute coordinates and first
    #         transforms them in 4D world coordinates, then takes the unit vector
    #         of that atom position to generate the loop animation. To generate
    #         the animation, first obtains the distance from the zero reference
    #         point (always 0,0,0) to the atom, then divides this distance in a
    #         defined number of cycles, this result will be the step for
    #         translation. For the translation, the world will move a number of
    #         steps defined, and every new point will be finded by multiplying the
    #         unit vector by the step. As a final step, to avoid biases, the world
    #         will be translated to the atom position in world coordinates.
    #         The effects will be applied on the model matrices of every VisMol
    #         object and the model matrix of the window.
    #     """
    #     if (self.zero_reference_point[0] != target[0]) or \
    #        (self.zero_reference_point[1] != target[1]) or \
    #        (self.zero_reference_point[2] != target[2]):
    #         self.zero_reference_point[:] = target
    #         pos = np.array([target[0],target[1],target[2],1], dtype=np.float32)
    #         model_pos = vismol_object.model_mat.T.dot(pos)[:3]
    #         self.model_mat = mop.my_glTranslatef(self.model_mat, -model_pos)
    #         unit_vec = model_pos / np.linalg.norm(model_pos)
    #         step = np.linalg.norm(model_pos)/15.0
    #         for i in range(15):
    #             to_move = unit_vec * step
                
    #             for index, vm_object in self.vm_session.vm_objects_dic.items():
    #                 vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -to_move)
                
    #             # WARNING: Method only works with GTK!!!
    #             if self.vm_session.toolkit == "Gtk_3.0":
    #                 self.parent_widget.get_window().invalidate_rect(None, False)
    #                 self.parent_widget.get_window().process_updates(False)
    #             elif self.vm_session.toolkit == "Qt5":
    #                 logger.critical("Not implemented for Qt5 yet :(")
    #                 raise RuntimeError("Not implemented for Qt5 yet :(")
    #             # WARNING: Method only works with GTK!!!
    #             time.sleep(self.vm_config.gl_parameters["center_on_coord_sleep_time"])
            
    #         for index, vm_object in self.vm_session.vm_objects_dic.items():
    #             model_pos = vm_object.model_mat.T.dot(pos)[:3]
    #             vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -model_pos)
            
    #         self.parent_widget.queue_draw()
    
