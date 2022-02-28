#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  representations.py
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

import ctypes
import numpy as np
from OpenGL import GL
from logging import getLogger

logger = getLogger(__name__)


class Representation:
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, name, active, indexes):
        self.vm_object = vismol_object
        self.vm_session = vismol_object.vm_session
        self.vm_glcore = vismol_glcore
        self.name = name
        self.active = active
        self.indexes = np.array(indexes, dtype=np.uint32)
        
        # representation
        self.vao = None
        self.ind_vbo = None
        self.coord_vbo = None
        self.col_vbo = None
        self.size_vbo = None
        # selection
        self.sel_vao = None
        self.sel_ind_vbo = None
        self.sel_coord_vbo = None
        self.sel_col_vbo = None
        self.sel_size_vbo = None
        # shaders
        self.shader_program = None
        self.sel_shader_program = None
    
    def _make_gl_vao_and_vbos(self):
        """ Function doc """
        self.shader_program = self.vm_glcore.shader_programs[self.name]
        self.sel_shader_program = self.vm_glcore.shader_programs[self.name + "_sel"]
        self._make_gl_representation_vao_and_vbos()
        self._make_gl_sel_representation_vao_and_vbos()
    
    def _make_gl_representation_vao_and_vbos(self, dot_sizes=None, normals=None):
        """ Function doc """
        self.vao       = self._make_gl_vao()
        self.ind_vbo   = self._make_gl_index_buffer(self.indexes)
        self.coord_vbo = self._make_gl_coord_buffer(self.vm_object.frames[0], self.shader_program)
        self.col_vbo   = self._make_gl_color_buffer(self.vm_object.colors, self.shader_program)
        if dot_sizes is not None:
            self.col_vbo = self._make_gl_color_buffer(colors, self.shader_program)
        if normals is not None and self.name == "surface":
            self.norm_vbo = self._make_gl_normal_buffer(normals, self.shader_program)
    
    def _make_gl_sel_representation_vao_and_vbos(self, dot_sizes=None, normals=None):
        """ Function doc """
        logger.debug("building '{}' background selection VAO and VBOs".format(self.name))
        self.sel_vao       = self._make_gl_vao()
        self.sel_ind_vbo   = self._make_gl_index_buffer(self.indexes)
        self.sel_coord_vbo = self._make_gl_coord_buffer(self.vm_object.frames[0], self.sel_shader_program)
        self.sel_col_vbo   = self._make_gl_color_buffer(self.vm_object.color_indexes, self.sel_shader_program)
        if dot_sizes is not None:
            self.sel_size_vbo = self._make_gl_size_buffer(dot_sizes ,self.sel_shader_program)
    
    def _make_gl_vao(self):
        """ Function doc """
        vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vao)
        return vao
    
    def _make_gl_index_buffer(self, indexes):
        """ Function doc """
        ind_vbo = GL.glGenBuffers(1)
        self.vm_session.vm_vbos.append(ind_vbo)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.nbytes, indexes, GL.GL_DYNAMIC_DRAW)
        return ind_vbo
    
    def _make_gl_coord_buffer(self, coords, program):
        """ Function doc """
        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_STATIC_DRAW)
        att_position = GL.glGetAttribLocation(program, "vert_coord")
        GL.glEnableVertexAttribArray(att_position)
        GL.glVertexAttribPointer(att_position, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        return coord_vbo
    
    def _make_gl_color_buffer(self, colors, program):
        """ Function doc """
        col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        att_colors = GL.glGetAttribLocation(program, "vert_color")
        GL.glEnableVertexAttribArray(att_colors)
        GL.glVertexAttribPointer(att_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
        return col_vbo
    
    def _make_gl_normal_buffer(self, normals, program):
        """ Function doc """
        normal_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, normal_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, normals.nbytes, normals, GL.GL_STATIC_DRAW)
        att_normals = GL.glGetAttribLocation(program, "vert_normal")
        if att_normals > 0:
            GL.glEnableVertexAttribArray(att_normals)
            GL.glVertexAttribPointer(att_normals, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*normals.itemsize, ctypes.c_void_p(0))
        return normal_vbo
    
    def _make_gl_size_buffer(self, dot_sizes, program):
        """ Function doc """
        size_vbo = GL.glGenBuffers(1)
        self.vm_session.vm_vbos.append(size_vbo)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, size_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, dot_sizes.nbytes, dot_sizes, GL.GL_STATIC_DRAW)
        att_size = GL.glGetAttribLocation(program, "vert_dot_size")
        GL.glEnableVertexAttribArray(att_size)
        GL.glVertexAttribPointer(att_size, 1, GL.GL_FLOAT, GL.GL_FALSE, dot_sizes.itemsize, ctypes.c_void_p(0))
        return size_vbo
    
    def _set_colors_to_buffer(self, col_vbo=True):
        """ Function doc """
        frame = self.vm_object.colors
        if col_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes, frame, GL.GL_STATIC_DRAW)
    
    def _set_coordinates_to_buffer(self, coord_vbo=True, sel_coord_vbo=True):
        """ This function assigns the coordinates to 
        be drawn by the function  draw_representation"""
        frame = self.vm_glcore._safe_frame_coords(self.vm_object)
        if coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes, frame, GL.GL_STATIC_DRAW)
        
        if sel_coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_coord_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes, frame, GL.GL_STATIC_DRAW)
    
    def _check_vao_and_vbos(self, indexes=None):
        """ Function doc """
        if self.sel_vao is None:
            self._make_gl_vao_and_vbos()
    
    def _enable_anti_alias_to_lines(self):
        """ Function doc """
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glEnable(GL.GL_LINE_SMOOTH)
        GL.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST)
    
    def _disable_anti_alias_to_lines(self):
        """ Function doc """
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
    
    def define_new_indexes_to_vbo(self, input_indexes=None):
        """ Function doc """
        if input_indexes is None:
            input_indexes = []
        self._check_vao_and_vbos()
        self.indexes = np.array(input_indexes, dtype=np.uint32)
        
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
        
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.sel_ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
    
    def change_vbo_colors(self, colors=None):
        """ Function doc """
        if colors is None:
            colors = []
        colors = np.array(colors, dtype=np.float32)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        att_colors = GL.glGetAttribLocation(self.shader_program, "vert_color")
        GL.glEnableVertexAttribArray(att_colors)
        GL.glVertexAttribPointer(att_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))


class DotsRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, indexes, active=True):
        """ Class initialiser """
        super(DotsRepresentation, self).__init__(vismol_object, vismol_glcore, "dots", active, indexes)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        height = self.vm_glcore.height
        GL.glUseProgram(self.shader_program)
        GL.glPointSize(1.5*height/abs(self.vm_glcore.dist_cam_zrp))
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.vm_object.atoms)), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glPointSize(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(self.sel_shader_program)
        GL.glPointSize(200/abs(self.vm_glcore.dist_cam_zrp))
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=False, sel_coord_vbo=True)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.vm_object.atoms)), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glPointSize(1)
        GL.glUseProgram(0)


class LinesRepresentation(Representation):
    """ Class doc """
    
    def __init__(self, vismol_object, vismol_glcore, indexes, active=True):
        """ Class initialiser """
        super(LinesRepresentation, self).__init__(vismol_object, vismol_glcore, "lines", active, indexes)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        line_width = self.vm_session.vm_config.gl_parameters["line_width"]
        line_width = (line_width*200/abs(self.vm_glcore.dist_cam_zrp)/2)**0.5
        GL.glLineWidth(line_width)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
            GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)

    def draw_background_sel_representation(self, line_width_factor=5):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        line_width = self.vm_session.vm_config.gl_parameters["line_width_selection"] 
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width) # line_width_factor -> turn the lines thicker
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=False, sel_coord_vbo=True)
            GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)


class NonBondedRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, indexes, active=True):
        """ Class initialiser """
        super(NonBondedRepresentation, self).__init__(vismol_object, vismol_glcore, "nonbonded", active, indexes)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        line_width = self.vm_session.vm_config.gl_parameters["line_width"]
        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(line_width*20/abs(self.vm_glcore.dist_cam_zrp))
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.vm_object.non_bonded_atoms)), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self, line_width_factor=5):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(20)
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=False, sel_coord_vbo=True)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.vm_object.non_bonded_atoms)), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)



'''
class DynamicBonds(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, name="sticks", indexes=None, active=True):
        """ Class initialiser """
        super(DynamicBonds, self).__init__(vismol_object, vismol_glcore, name, active, indexes)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos ()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(40/abs(self.vm_glcore.dist_cam_zrp))
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        self.vm_glcore.load_lights(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            frame = self.vm_glcore.frame
            if frame < len(self.vm_object.dynamic_bons):
                self.define_new_indexes_to_vbo(self.vm_object.dynamic_bons[frame])
                self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.dynamic_bons[frame])*2), GL.GL_UNSIGNED_INT, None)
            else:
                self.define_new_indexes_to_vbo(self.vm_object.dynamic_bons[-1])
                self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.dynamic_bons[-1])*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(20)
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            frame = self.vm_glcore.frame
            if frame < len(self.vm_object.dynamic_bons):
                self.define_new_indexes_to_vbo(self.vm_object.dynamic_bons[frame])
                self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.dynamic_bons[frame])*2), GL.GL_UNSIGNED_INT, None)
            else:
                self.define_new_indexes_to_vbo(self.vm_object.dynamic_bons[-1])
                self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.dynamic_bons[-1])*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)


class SticksRepresentation(Representation):
    """ Class doc """
    
    def __init__(self, vismol_object, vismol_glcore, name="sticks", indexes=None, active=True):
        """ Class initialiser """
        super(SticksRepresentation, self).__init__(vismol_object, vismol_glcore, name, active, indexes)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos ()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(40/abs(self.vm_glcore.dist_cam_zrp))
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        self.vm_glcore.load_lights(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
            GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(0)
        GL.glLineWidth(1)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(20)
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=False, sel_coord_vbo=True)
            GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)


class RibbonsRepresentation(Representation):
    """ Class doc """
    
    def __init__(self, vismol_object, vismol_glcore, name="ribbon", indexes=None, active=True):
        """ Class initialiser """
        super(RibbonsRepresentation, self).__init__(vismol_object, vismol_glcore, name, active, indexes)
        
        if self.vm_object.c_alpha_bonds == []:
            self.vm_object.get_backbone_indexes()
        
        indexes = []
        for bond in self.vm_object.c_alpha_bonds:
            indexes.append(bond.atom_index_i)
            indexes.append(bond.atom_index_j)
        
        if indexes == []:
            self.activate = False
        else:
            self.indexes = np.array(indexes, dtype=np.uint32)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos ()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        ribbon_width = self.vm_session.vm_config.gl_parameters["ribbon_width"]
        ribbon_width = (ribbon_width*20)/abs(self.vm_glcore.dist_cam_zrp)/2
        GL.glLineWidth(ribbon_width)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
            GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        line_width = self.vm_session.vm_config.gl_parameters["line_width_selection"] 
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width)
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=False, sel_coord_vbo=True)
            GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)


class SpheresRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, name="dots", indexes=None, active=True):
        """ Class initialiser """
        super(SpheresRepresentation, self).__init__(vismol_object, vismol_glcore, name, active, indexes)
        
        self.level= self.vm_session.vm_config.gl_parameters['sphere_quality']
        self.scale= self.vm_session.vm_config.gl_parameters['sphere_scale']
    
    def update_atomic_indexes(self, indexes=None):
        """ Function doc """
        for index in indexes:
            self.atoms.append(self.vm_object.atoms[index])
            self.atomic_indexes.append(index)
        self._create_sphere_data()
        self._update_sphere_data_to_vbos()
        self.active = True
        
    def _update_sphere_data_to_vbos(self):
        """ Function doc """
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
        
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.nbytes, self.coords, GL.GL_STATIC_DRAW)
        
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.centr_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
        
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.colors.itemsize*len(self.colors), self.colors, GL.GL_STATIC_DRAW)
    
    def _create_sphere_data(self):
        """ Function doc """
        init = time.time()
        qtty = len(self.atoms)
        nucleus = [0.0, 0.0, 0.0] * qtty
        colores = [0.0, 0.0, 0.0] * qtty
        coords  = sphd.sphere_vertices[self.level] * qtty
        centers = sphd.sphere_vertices[self.level] * qtty
        colors  = sphd.sphere_vertices[self.level] * qtty
        indexes = np.array(sphd.sphere_triangles[self.level] * qtty, dtype=np.uint32)
        elems  = int(len(sphd.sphere_vertices[self.level]) / 3)
        offset = int(len(sphd.sphere_vertices[self.level]))
        inds_e = int(len(sphd.sphere_triangles[self.level]))
        
        self.centers_list = []
        self.frames = []
        frame =0
        
        for a, atom in enumerate(self.atoms):
            pos = atom.coords(frame)
            colors[a*offset:(a+1)*offset]  = [atom.color[0],atom.color[1],atom.color[2]]*elems
            centers[a*offset:(a+1)*offset] = [pos[0],pos[1],pos[2]]*elems
            for i in range(elems):
                coords[a*offset+i*3]   *= atom.radius * self.scale
                coords[a*offset+i*3+1] *= atom.radius * self.scale
                coords[a*offset+i*3+2] *= atom.radius * self.scale
                coords[a*offset+i*3]   += pos[0]
                coords[a*offset+i*3+1] += pos[1]
                coords[a*offset+i*3+2] += pos[2]
            indexes[a*inds_e:(a+1)*inds_e] += a*elems
        end = time.time()
        print('Time used creating nucleus, vertices and colors:', end-init)
        
        self.coords = np.array(coords, dtype=np.float32)
        self.frames.append(self.coords)
        self.centers = np.array(centers, dtype=np.float32)
        self.centers_list.append(self.centers)
        self.colors  = np.array(colors, dtype=np.float32)
        self.indexes = indexes
        self.triangles = len(self.indexes)
        init = time.time()
        if len(self.vm_object.frames) > 1:
            for frame in range(1, len(self.vm_object.frames) - 1):
                coords  = sphd.sphere_vertices[self.level]*qtty
                centers = sphd.sphere_vertices[self.level]*qtty
                for a, atom in enumerate(self.atoms):
                    pos = atom.coords(frame)
                    centers[a*offset:(a+1)*offset] = [pos[0],pos[1],pos[2]]*elems
                    for i in range(elems):
                        coords[a*offset+i*3]   *= atom.radius * self.scale
                        coords[a*offset+i*3+1] *= atom.radius * self.scale
                        coords[a*offset+i*3+2] *= atom.radius * self.scale
                        coords[a*offset+i*3]   += pos[0]
                        coords[a*offset+i*3+1] += pos[1]
                        coords[a*offset+i*3+2] += pos[2]
                self.coords  = np.array(coords, dtype=np.float32)
                self.frames.append(self.coords)
                self.centers = np.array(centers, dtype=np.float32)
                self.centers_list.append(self.centers)
        end = time.time()
        print('Time used creating nucleus, vertices and colors:', end-init)
    
    def _make_gl_vao_and_vbos(self):
        """ Function doc """
        self.shader_program     = self.vm_glcore.shader_programs[self.name]
        self.sel_shader_program = self.vm_glcore.shader_programs[self.name + '_sel']
        
        self.vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(self.vao)
        self.ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
        
        self.coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.nbytes, self.coords, GL.GL_STATIC_DRAW)
        gl_coord = GL.glGetAttribLocation(self.shader_program , 'vert_coord')
        GL.glEnableVertexAttribArray(gl_coord)
        GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.coords.itemsize, ctypes.c_void_p(0))
        
        self.centr_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.centr_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
        gl_center = GL.glGetAttribLocation(self.shader_program , 'vert_centr')
        GL.glEnableVertexAttribArray(gl_center)
        GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.centers.itemsize, ctypes.c_void_p(0))
        
        self.col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.colors.itemsize*len(self.colors), self.colors, GL.GL_STATIC_DRAW)
        gl_colors = GL.glGetAttribLocation(self.shader_program, 'vert_color')
        GL.glEnableVertexAttribArray(gl_colors)
        GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.colors.itemsize, ctypes.c_void_p(0))
        
        self.triangles = len(self.indexes)
        
        colors_idx = self.vm_object.color_indexes
        self._make_gl_sel_representation_vao_and_vbos(indexes=self.indexes, coords=self.coords,
                                                      colors=colors_idx, dot_sizes=None)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        GL.glUseProgram(self.shader_program)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_lights(self.shader_program)
        self.vm_glcore.load_fog(self.shader_program)
        
        if self.vao is not None:
            GL.glBindVertexArray (self.vao)
            if self.vm_glcore.modified_view:
                pass
            else:
                frame = self.vm_glcore._get_vismol_object_frame (self.vm_object)
                self.coords = self.frames[frame]
                self.centers =self.centers_list[frame]
                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
                GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.itemsize*int(len(self.coords)), self.coords, GL.GL_STATIC_DRAW)
                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.centr_vbo)
                GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
                GL.glDrawElements(GL.GL_TRIANGLES,  self.triangles , GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
    
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass


class ImpostorRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, name = "impostor", indexes=None, active=True, scale=1.0):
        """ Class initialiser """
        super(ImpostorRepresentation, self).__init__(vismol_object, vismol_glcore, name, active, indexes)
        self.scale = scale
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        height = self.vm_glcore.height
        dist_cam_zrp = self.vm_glcore.dist_cam_zrp
        xyz_coords = self.vm_glcore.glcamera.get_modelview_position(self.vm_object.model_mat)
        u_campos = GL.glGetUniformLocation(self.shader_program, "u_campos")
        GL.glUniform3fv(u_campos, 1, xyz_coords)
        self.vm_glcore.load_lights(self.shader_program)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
            GL.glDrawElements(GL.GL_POINTS, len(self.vm_object.atoms), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glPointSize(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(self.sel_shader_program)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=False, sel_coord_vbo=True)
            GL.glDrawElements(GL.GL_POINTS, len(self.vm_object.atoms), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glPointSize(1)
        GL.glUseProgram(0)


class CartoonRepresentation(Representation):
    def __init__ (self, name = 'cartoon', active = True, rep_type = 'mol', vismol_object = None, vismol_glcore = None, indexes = []):
        self.name               = name
        self.active             = active
        self.type               = rep_type

        self.vm_object             = vismol_object
        self.vm_glcore             = vm_glcore
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.norm_vbo       = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        
        
        coords, normals, indexes, colors = cartoon.cartoon(vismol_object, spline_detail=5)
        
        coords = coords.flatten()
        normals = normals.flatten()
        colors = colors.flatten()
        
        
        self.coords2 = coords
        self.colors2 = colors
        self.normals2 = normals
        self.indexes2 = indexes


    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        #if indexes is not None:
        #    pass
        #else:
        
        #dot_qtty  = int(len(self.vm_object.frames[0])/3)
        #indexes = []
        #for i in range(dot_qtty):
        #    indexes.append(i)
        

        self.shader_program     = self.vm_glcore.shader_programs[self.name]
        #self.sel_shader_program = self.vm_glcore.shader_programs[self.name+'_sel']
        

        """
        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.array(self.colors2, dtype=np.float32)
        normals = np.array(self.normals2, dtype=np.float32)
        indexes = np.array(self.indexes2, dtype=np.uint32)
        """
        
        
        coords  = self.coords2 
        colors  = self.colors2 
        normals = self.normals2
        indexes = self.indexes2
        
        print ('len(coords),len(colors), len(normals),len(indexes)', len(coords),len(colors), len(normals),len(indexes)  )

        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   normals    = normals
                                                   )
        
        
        
        self.ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.itemsize*len(indexes), indexes, GL.GL_DYNAMIC_DRAW)
        
        #self.coord_vbo = GL.glGenBuffers(1)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
        ##GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
        #GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_STATIC_DRAW)
        #gl_coord = GL.glGetAttribLocation(self.shader_program, 'vert_coord')
        #GL.glEnableVertexAttribArray(gl_coord)
        #GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        
        
        self.col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.itemsize*len(colors), colors, GL.GL_STATIC_DRAW)
        gl_color = GL.glGetAttribLocation(self.shader_program, 'vert_color')
        GL.glEnableVertexAttribArray(gl_color)
        GL.glVertexAttribPointer(gl_color, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))

        self.norm_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.norm_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, normals.itemsize*len(normals), normals, GL.GL_STATIC_DRAW)
        gl_norm = GL.glGetAttribLocation(self.shader_program, 'vert_norm')
        GL.glEnableVertexAttribArray(gl_norm)
        GL.glVertexAttribPointer(gl_norm, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*normals.itemsize, ctypes.c_void_p(0))
        
        
        
        
        
        #self.centr_vbo = GL.glGenBuffers(1)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coords)
        #GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
        #gl_center = GL.glGetAttribLocation(self.shader_program , 'vert_centr')
        #GL.glEnableVertexAttribArray(gl_center)
        #GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        
        
        
        colors_idx = self.vm_object.color_indexes
        self.sel_vao = True
        """
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )
        """
    def draw_representation (self):
        """ Function doc """
        self._check_vao_and_vbos ()
        #self._enable_anti_alias_to_lines()
        
        
        
        
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glDisable(GL.GL_CULL_FACE)
        #GL.glCullFace(GL.GL_BACK)
        view = self.vm_glcore.glcamera.view_matrix
        
        GL.glUseProgram(self.shader_program )
        
        #print (self.vm_object.model_mat,view)
        
        m_normal = np.array(np.matrix(np.dot(view, self.vm_object.model_mat)).I.T)
        
        self.vm_glcore.load_matrices(self.shader_program , self.vm_object.model_mat)
        self.vm_glcore.load_lights  (self.shader_program )
        self.vm_glcore.load_fog     (self.shader_program )
        GL.glBindVertexArray(self.vao)
        
        
        
        
        
        
        
        
        
        """
        #print ("DotsRepresentation")
        height = self.vm_glcore.height
        
        GL.glUseProgram(self.shader_program)
        #1*self.height dot_size
        #GL.glLineWidth(40/abs(self.vm_glcore.dist_cam_zrp))
        GL.glPointSize(0.1*height/abs(self.vm_glcore.dist_cam_zrp)) # dot size not included yet
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        """
        if self.vm_glcore.modified_view:
            pass
        
        else:
            """
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea"""
            # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
            #GL.glDrawElements(GL.GL_POINTS, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_LOOP, int(len(self.coords2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_STRIP, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            
            #print("int(len(self.indexes2))", int(len(self.indexes2)))
            GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_TRIANGLES, 54060, GL.GL_UNSIGNED_INT, None)
        
        #GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        #GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass


class SurfaceRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, name = "surface", active = True, rep_type = "mol", vismol_object = None, vm_glcore = None, indexes = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = rep_type

        self.vm_object             = vismol_object
        self.vm_glcore             = vm_glcore
        
        

        
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.norm_vbo       = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        self.read_surface_data()
    
    
    ##### sub 2 vev3 vectors
    def sub_vec3(self, a, b):
        c = [ a[0] - b[0],
              a[1] - b[1],
              a[2] - b[2] ]

        return c

    ## add 2 vectors and take the avg
    ## if a vector is still 0 we just take b
    def avg_add_vec3(self, a, b):
        if a[0] == 0.0 and a[1] == 0.0 and a[2] == 0.0 :
            return b

        c = [ (a[0] + b[0]) * 0.5 ,
              (a[1] + b[1]) * 0.5 ,
              (a[2] + b[2]) * 0.5 ]

        return c    

    ## make the cross product of 2 vectors
    def cross_vec3(self, a, b):
        c = [a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0]]

        return c
    #############################################
        
    
    
    def read_surface_data(self):
        """ Function doc """
        #from random import random 
        #
        #[verts, tris, verts_gpu, tris_gpu] = edtsurf.calc_surface("/home/fernando/programs/EasyHybrid3/Coords/pdbs/1bx4_H.pdb")
        #self.coords2  = verts_gpu
        #self.indexes2 = tris_gpu
        #self.colors2  = []
        #
        #
        #size = len( self.coords2 )
        #for i in range(size):
        #    self.colors2.append(float(i/size) + random())
        
        rawdata = open("../EasyHybrid3/Coords/pdbs/1bx4.ply", "r")
        lines  = rawdata.readlines()
        
        self.coords2 = []
        self.colors2 = []
        self.normals2 = []
        self.indexes2 = []
        avg_normals_indexes = []
        
        
        for line in lines:
            line2 = line.split()
            
            if len(line2) == 6:
                #print (line2)
                self.coords2.append(float(line2[0]))
                self.coords2.append(float(line2[1]))
                self.coords2.append(float(line2[2]))
                                                  
                self.colors2.append(float(line2[3])/255)
                self.colors2.append(float(line2[4])/255)
                self.colors2.append(float(line2[5])/255)
                
                self.normals2.append(float(line2[0]))
                self.normals2.append(float(line2[1]))
                self.normals2.append(float(line2[2]))                
                avg_normals_indexes.append( ( 0.0 , 0.0 , 0.0 ) )  ### NEW !!! 

            if len(line2) == 7:
                
                self.indexes2.append(int(line2[1]))
                self.indexes2.append(int(line2[2]))
                self.indexes2.append(int(line2[3]))
                
        
        ## calculate normals and interpolate them (thanks a lot Kai)
        for i in range( 0 , len(self.indexes2) , 3 ):

            index_1 = self.indexes2[i] * 3;
            index_2 = self.indexes2[i+1] * 3;
            index_3 = self.indexes2[i+2] * 3;
            vertex_1 = ( self.coords2[index_1] , self.coords2[index_1+1] , self.coords2[index_1+2] )
            vertex_2 = ( self.coords2[index_2] , self.coords2[index_2+1] , self.coords2[index_2+2] )
            vertex_3 = ( self.coords2[index_3] , self.coords2[index_3+1] , self.coords2[index_3+2] )

            vec_p0_p1 = self.sub_vec3( vertex_2 , vertex_1 )
            vec_p0_p2 = self.sub_vec3( vertex_3 , vertex_1 )
            norm_vec  = self.cross_vec3( vec_p0_p1, vec_p0_p2 )

            vert_index_1 = self.indexes2[i] ;
            vert_index_2 = self.indexes2[i+1] ;
            vert_index_3 = self.indexes2[i+2] ;
            
            avg_normals_indexes[vert_index_1] = self.avg_add_vec3( avg_normals_indexes[vert_index_1] , norm_vec )
            avg_normals_indexes[vert_index_2] = self.avg_add_vec3( avg_normals_indexes[vert_index_2] , norm_vec )
            avg_normals_indexes[vert_index_3] = self.avg_add_vec3( avg_normals_indexes[vert_index_3] , norm_vec )


        ## set all new interpolated normals   
        for i in range( 0 , len(self.indexes2) , 1 ):
            index_1 = self.indexes2[i] * 3;

            self.normals2[index_1]   = avg_normals_indexes[self.indexes2[i]][0]
            self.normals2[index_1+1] = avg_normals_indexes[self.indexes2[i]][1]
            self.normals2[index_1+2] = avg_normals_indexes[self.indexes2[i]][2]





               
                

    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        #if indexes is not None:
        #    pass
        #else:
        
        #dot_qtty  = int(len(self.vm_object.frames[0])/3)
        #indexes = []
        #for i in range(dot_qtty):
        #    indexes.append(i)
        

        self.shader_program     = self.vm_glcore.shader_programs[self.name]
        self.sel_shader_program = self.vm_glcore.shader_programs[self.name+"_sel"]
        
        #indexes = np.array(self.vm_object.index_bonds, dtype=np.uint32)
        #indexes = np.array(self.vm_object.idex, dtype=np.uint32)

        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.array(self.colors2, dtype=np.float32)
        normals = np.array(self.normals2, dtype=np.float32)
        #indexes = range(0, len(self.coords2))     
        #indexes = np.array(indexes, dtype=np.uint32)
        indexes = np.array(self.indexes2, dtype=np.uint32)
        #print (indexes)


        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   normals    = normals
                                                   )
        
        #self.centr_vbo = GL.glGenBuffers(1)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coords)
        #GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
        #gl_center = GL.glGetAttribLocation(self.shader_program , "vert_centr")
        #GL.glEnableVertexAttribArray(gl_center)
        #GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        
        
        
        colors_idx = self.vm_object.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_vao_and_vbos ()
        #self._enable_anti_alias_to_lines()
        
        
        
        
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        view = self.vm_glcore.glcamera.view_matrix
        
        GL.glUseProgram(self.shader_program )
        
        #print (self.vm_object.model_mat,view)
        
        m_normal = np.array(np.matrix(np.dot(view, self.vm_object.model_mat)).I.T)
        
        self.vm_glcore.load_matrices(self.shader_program , self.vm_object.model_mat)
        self.vm_glcore.load_lights  (self.shader_program )
        self.vm_glcore.load_fog     (self.shader_program )
        GL.glBindVertexArray(self.vao)
        
        
        
        
        
        
        
        
        
        """
        #print ("DotsRepresentation")
        height = self.vm_glcore.height
        
        GL.glUseProgram(self.shader_program)
        #1*self.height dot_size
        #GL.glLineWidth(40/abs(self.vm_glcore.dist_cam_zrp))
        GL.glPointSize(0.1*height/abs(self.vm_glcore.dist_cam_zrp)) # dot size not included yet
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        """
        if self.vm_glcore.modified_view:
            pass
        
        else:
            """
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea"""
            # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
            #GL.glDrawElements(GL.GL_POINTS, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_LOOP, int(len(self.coords2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_STRIP, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
        
        #GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        #GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass


class WiresRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, name = "wires", active = True, rep_type = "mol", vismol_object = None, vm_glcore = None, indexes = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = rep_type
        self.vm_object             = vismol_object
        self.vm_glcore             = vm_glcore
        
        # representation    
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
        
        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None

        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        self.read_surface_data()
    
    def read_surface_data(self):
        """ Function doc """
        rawdata = open("../EasyHybrid3/Coords/pdbs/1bx4.ply", "r")
        lines  = rawdata.readlines()
        
        self.coords2 = []
        self.colors2 = []
        self.indexes2 = []
        
        for line in lines:
            line2 = line.split()
            if len(line2) == 6:
                self.coords2.append(float(line2[0]))
                self.coords2.append(float(line2[1]))
                self.coords2.append(float(line2[2]))
                self.colors2.append(float(line2[3])/255)
                self.colors2.append(float(line2[4])/255)
                self.colors2.append(float(line2[5])/255)
            if len(line2) == 7:
                self.indexes2.append(int(line2[1]))
                self.indexes2.append(int(line2[2]))
                self.indexes2.append(int(line2[3]))

    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        self.shader_program     = self.vm_glcore.shader_programs[self.name]
        self.sel_shader_program = self.vm_glcore.shader_programs[self.name+"_sel"]
        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.zeros(len(self.colors2))
        indexes = np.array(self.indexes2, dtype=np.uint32)
        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.vm_object.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_vao_and_vbos ()
        pass
        #GL.glEnable(GL.GL_DEPTH_TEST)
        #GL.glEnable(GL.GL_CULL_FACE)
        #GL.glCullFace(GL.GL_BACK)
        #
        ##LineWidth = (80/abs(self.vm_glcore.dist_cam_zrp)/2)**0.5  #40/abs(self.vm_glcore.dist_cam_zrp)
        ##GL.glLineWidth(2)
        #
        #
        #view = self.vm_glcore.glcamera.view_matrix
        #GL.glUseProgram(self.shader_program )
        #m_normal = np.array(np.matrix(np.dot(view, self.vm_object.model_mat)).I.T)
        #self.vm_glcore.load_matrices(self.shader_program , self.vm_object.model_mat)
        ##self.vm_glcore.load_lights  (self.shader_program )
        #self.vm_glcore.load_fog     (self.shader_program )
        #GL.glBindVertexArray(self.vao)
        #if self.vm_glcore.modified_view:
        #    pass
        #
        #else:
        #    """
        #    This function checks if the number of the called frame will not exceed 
        #    the limit of frames that each object has. Allowing two objects with 
        #    different trajectory sizes to be manipulated at the same time within the 
        #    glArea"""
        #    # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
        #    GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
        #GL.glDisable(GL.GL_DEPTH_TEST)
        
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass


class LabelRepresentation:
    """ Class doc """
    
    def __init__ (self, name = "labels", active = True, rep_type = "mol", vismol_object = None, vm_glcore = None, indexes = []):
        """ Class initialiser """
        self.vm_object = vismol_object
        self.name   = name
        self.active = True
        self.vm_glcore = vm_glcore
        
        self.chars     = 0 
        #self._check_vao_and_vbos()
        
    def _check_vao_and_vbos (self, indexes = None):
        """ Function doc """
        if self.vm_object.vm_font.vao is None:
            self.vm_object.vm_font.make_freetype_font()
            self.vm_object.vm_font.make_freetype_texture(self.vm_glcore.freetype_program)
        
        if self.chars == 0:
            print("self._build_buffer()")
            self._build_buffer()
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vm_object.vm_font.vbos[0])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.xyz_pos.itemsize*len(self.xyz_pos), self.xyz_pos, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vm_object.vm_font.vbos[1])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.uv_coords.itemsize*len(self.uv_coords), self.uv_coords, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)


    def _build_buffer (self, indexes = None):
        self.chars     = 0
        self.xyz_pos   = []
        self.uv_coords = []
        for atom in self.vm_object.atoms:
            
            texto = atom.name
            point = np.array(atom.coords (self.vm_glcore.frame),np.float32)
            point = np.array((point[0],point[1],point[2],1),np.float32)
            point = np.dot(point, self.vm_object.model_mat)

            GL.glBindTexture(GL.GL_TEXTURE_2D, self.vm_object.vm_font.texture_id)
            for i,c in enumerate(texto):
                self.chars += 1
                c_id = ord(c)
                x = c_id%16
                y = c_id//16-2
                self.xyz_pos.append(point[0]+i*self.vm_object.vm_font.char_width)
                self.xyz_pos.append(point[1])
                self.xyz_pos.append(point[2])

                self.uv_coords.append(x*self.vm_object.vm_font.text_u)
                self.uv_coords.append(y*self.vm_object.vm_font.text_v)
                self.uv_coords.append((x+1)*self.vm_object.vm_font.text_u)
                self.uv_coords.append((y+1)*self.vm_object.vm_font.text_v)
            #print(texto)
        #print("xyz_pos  ",len(self.xyz_pos))
        #print("uv_coords",len(self.uv_coords))
        #print("atoms    ",len(self.vm_object.atoms))
        #print("chars    ",self.chars)
        
        self.xyz_pos   = np.array(self.xyz_pos  , np.float32)
        self.uv_coords = np.array(self.uv_coords, np.float32)
        


    
    
    def draw_representation (self):
        """ Function doc """
        self._check_vao_and_vbos()
        
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glUseProgram(self.vm_glcore.freetype_program)
        
        self.vm_object.vm_font.load_matrices(self.vm_glcore.freetype_program, self.vm_glcore.glcamera.view_matrix, self.vm_glcore.glcamera.projection_matrix)
        self.vm_object.vm_font.load_font_params(self.vm_glcore.freetype_program)
        
        GL.glBindVertexArray(self.vm_object.vm_font.vao)
        GL.glDrawArrays(GL.GL_POINTS, 0, self.chars)
        GL.glDisable(GL.GL_BLEND)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)

        

    def draw_background_sel_representation  (self):
        """ Function doc """
        pass

'''