#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from OpenGL import GL
from logging import getLogger
from vismol.libgl.glaxis import GLAxis
from vismol.libgl.glcamera import GLCamera
import vismol.utils.matrix_operations as mop


logger = getLogger(__name__)


class VismolGLCore():
    
    def __init__(self, widget, vismol_config, width, height):
        self.parent_widget = widget
        self.vm_config = vismol_config
        self.shader_programs = {}
        self.width = np.float32(width)
        self.height = np.float32(height)
    
    def initialize_parent(self):
        """ Enables the buffers and other charasteristics of the OpenGL context.
            sets the initial projection, view and model matrices.
        """
        self.model_mat = np.identity(4, dtype=np.float32)
        self.zero_reference_point = np.zeros(3, dtype=np.float32)
        self.glcamera = GLCamera(self.vm_config.gl_parameters["field_of_view"],
                                 self.width / self.height,
                                 np.array([0, 0, 10], dtype=np.float32),
                                 self.zero_reference_point)
        self.axis = GLAxis(vm_glcore = self)
        self.parent_widget.set_has_depth_buffer(True)
        self.parent_widget.set_has_alpha(True)
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
        self.mouse_button = None
        self.dist_cam_zrp = np.linalg.norm(self.glcamera.get_position())
        self.shader_flag = True
        self.dragging = False
        self.show_axis = True
        self.mouse_x, self.mouse_y = np.float32(0.0), np.float32(0.0)
        self.mouse_rotate, self.mouse_zoom, self.mouse_pan = False, False, False
        self.drag_pos_x, self.drag_pos_y, self.drag_pos_z = None, None, None
    
    def resize_window(self, width, height):
        """ Resizing function, takes the widht and height of the widget
            and modifies the view in the camera acording to the new values
        
            Keyword arguments:
            width -- Actual width of the window
            height -- Actual height of the window
        """
        self.width = np.float32(width)
        self.height = np.float32(height)
        self.right = self.width / self.height
        self.left = -self.right
        self.center_x = self.width / 2.0
        self.center_y = self.height / 2.0
        self.glcamera.viewport_aspect_ratio = self.width / self.height
        _proj_mat = mop.my_glPerspectivef(self.glcamera.field_of_view,
                                          self.glcamera.viewport_aspect_ratio,
                                          self.glcamera.z_near, self.glcamera.z_far)
        self.glcamera.set_projection_matrix(_proj_mat)
    
    def mouse_pressed(self, button_number, mouse_x, mouse_y):
        logger.critical("NotImplementedError, the child class must implement mouse_pressed")
        raise NotImplementedError("Subclasses must implement this method")
    
    def mouse_released(self, button_number, mouse_x, mouse_y):
        logger.critical("NotImplementedError, the child class must implement mouse_released")
        raise NotImplementedError("Subclasses must implement this method")
    
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
            self.parent_widget.queue_draw()
    
    def get_viewport_pos(self, x, y):
        """ Function doc """
        px = (2.0*x - self.width) / self.width
        py = (self.height - 2.0*y) / self.height
        return [px, py, self.glcamera.z_near]
    
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
    
    def queue_draw(self):
        """ Function doc """
        self.parent_widget.queue_draw()
    
    def mouse_scroll(self, direction):
        logger.critical("NotImplementedError, the child class must implement mouse_scroll")
        raise NotImplementedError("Subclasses must implement this method")
    
    def _rotate_view(self, dx, dy, x, y):
        logger.critical("NotImplementedError, the child class must implement _rotate_view")
        raise NotImplementedError("Subclasses must implement this method")
    
    def _pan_view(self, x, y):
        logger.critical("NotImplementedError, the child class must implement _pan_view")
        raise NotImplementedError("Subclasses must implement this method")
    
    def _zoom_view(self, dy):
        logger.critical("NotImplementedError, the child class must implement _zoom_view")
        raise NotImplementedError("Subclasses must implement this method")
    
    def render(self):
        logger.critical("NotImplementedError, the child class must implement render")
        raise NotImplementedError("Subclasses must implement this method")
    
    def create_gl_programs(self):
        logger.critical("NotImplementedError, the child class must implement create_gl_programs")
        raise NotImplementedError("Subclasses must implement this method")
    
    def _safe_frame_coords(self):
        logger.critical("NotImplementedError, the child class must implement _safe_frame_coords")
        raise NotImplementedError("Subclasses must implement this method")
    
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
            logger.critical("Error compiling the shader: {}".format(shader_type))
            raise RuntimeError(logger.critical(GL.glGetShaderInfoLog(shader)))
        return shader
    
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
        
        # it is not necessary get the location everytime - change it later pls!
        model = GL.glGetUniformLocation(program, "model_mat")
        GL.glUniformMatrix4fv(model, 1, GL.GL_FALSE, model_mat)
        view = GL.glGetUniformLocation(program, "view_mat")
        GL.glUniformMatrix4fv(view, 1, GL.GL_FALSE, self.glcamera.view_matrix)
        proj = GL.glGetUniformLocation(program, "proj_mat")
        GL.glUniformMatrix4fv(proj, 1, GL.GL_FALSE, self.glcamera.projection_matrix)
    
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

