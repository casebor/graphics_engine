#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import math
import numpy as np

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

import OpenGL
from OpenGL import GL

import libs.vaos as vaos
import libs.camera as cam

class VMWindow(Gtk.GLArea):
    """docstring for VMWindow"""
    def __init__(self, width=350, height=250, **kwargs):
        super().__init__()
        # self.gl_area = Gtk.GLArea.new()
        self.connect("realize", self.initialize)
        self.connect("render", self.render)
        self.connect("resize", self.reshape_window)
        self.connect("key-press-event", self.key_pressed)
        self.connect("key-release-event", self.key_released)
        self.connect("scroll-event", self.mouse_scroll)
        self.connect("button-press-event", self.mouse_pressed)
        self.connect("motion-notify-event", self.mouse_motion)
        self.connect("button-release-event", self.mouse_released)
        self.grab_focus()
        self.set_events( self.get_events() | Gdk.EventMask.SCROLL_MASK
                       | Gdk.EventMask.BUTTON_PRESS_MASK | Gdk.EventMask.BUTTON_RELEASE_MASK
                       | Gdk.EventMask.POINTER_MOTION_MASK | Gdk.EventMask.POINTER_MOTION_HINT_MASK
                       | Gdk.EventMask.KEY_PRESS_MASK | Gdk.EventMask.KEY_RELEASE_MASK)
        self.vert_shader = kwargs.get("vertex_shader", None)
        self.frag_shader = kwargs.get("fragment_shader", None)
        self.geom_shader = kwargs.get("geometry_shader", None)
        self.draw_type = kwargs.get("draw_type", "circles")
        self.gl_data = kwargs.get("data", None)
        self.gl_programs_compiled = False
        self.gl_program = None
        self.animation = False
    
    def initialize(self, area):
        """ Function doc """
        aloc = self.get_allocation()
        self.width = float(aloc.width)
        self.height = float(aloc.height)
        self.right = self.width / self.height
        self.left = -self.right
        self.z_near = 0.1
        self.z_far = 100
        self.min_znear = 0.1
        self.min_zfar = 1.1
        self.scroll = 0.3
        self.fov = 20.0 # Field Of View = fov
        self.var = self.width/self.height # Viewport Aspect Ratio
        self.model_mat = np.identity(4, dtype=float)
        self.view_mat = cam.my_glTranslatef(np.identity(4, dtype=float), [0, 0, -1])
        self.proj_mat = cam.my_glOrthof(self.width, self.height, self.z_near, self.z_far)
        self.bckgrnd_color = np.array([0, 0, 0, 1], dtype=float)
        self.vao = None
        self.elements = None
        self.ctrl = False
        self.cam_pos = self.get_cam_pos()
        self.mouse_x = 0
        self.mouse_y = 0
        self.mouse_rotate = False
        self.mouse_zoom = False
        self.mouse_pan = False
        self.light_position = np.array([-2.5, 2.5, 2.5],dtype=float)
        self.light_color = np.array([1.0, 1.0, 1.0, 1.0],dtype=float)
        self.light_ambient_coef = 0.5
        self.light_shininess = 5.5
        self.light_intensity = np.array([0.6, 0.6, 0.6],dtype=float)
        self.light_specular_color = np.array([1.0, 1.0, 1.0],dtype=float)
        self.make_current()
        if (self.get_error() != None):
            return
        self.set_has_depth_buffer(True)
        self.set_has_alpha(True)
        self._create_gl_programs()
    
    def get_cam_pos(self):
        """ Returns the position of the camera in XYZ coordinates
            The type of data returned is 'numpy.ndarray'.
        """
        modelview = cam.my_glMultiplyMatricesf(self.model_mat, self.view_mat)
        crd_xyz = -1 * np.mat(modelview[:3,:3]) * np.mat(modelview[3,:3]).T
        return crd_xyz.A1
    
    def _update_cam_pos(self):
        """ Function doc """
        self.cam_pos = self.get_cam_pos()
    
    def reshape_window(self, widget, width, height):
        """ Function doc """
        self.width = float(width)
        self.height = float(height)
        self.proj_mat = cam.my_glOrthof(self.width, self.height, self.z_near, self.z_far)
    
    def _create_gl_programs(self):
        """ Function doc """
        print('OpenGL version: ',GL.glGetString(GL.GL_VERSION))
        try:
            print('OpenGL major version: ',GL.glGetDoublev(GL.GL_MINOR_VERSION))
            print('OpenGL minor version: ',GL.glGetDoublev(GL.GL_MAJOR_VERSION))
        except:
            print('OpenGL major version not found')
        self.gl_program = self._load_shaders(self.vert_shader, self.frag_shader,
                                             self.geom_shader)
    
    def _load_shaders(self, vertex, fragment, geometry):
        """ Here the shaders are loaded and compiled to an OpenGL program. By default
            the constructor shaders will be used, if you want to change the shaders
            use this function. The flag is used to create only one OpenGL program.
            
            Keyword arguments:
            vertex -- The vertex shader to be used
            fragment -- The fragment shader to be used
        """
        my_vertex_shader = self._create_shader(vertex, GL.GL_VERTEX_SHADER)
        my_fragment_shader = self._create_shader(fragment, GL.GL_FRAGMENT_SHADER)
        if geometry is not None:
            my_geometry_shader = self._create_shader(geometry, GL.GL_GEOMETRY_SHADER)
        program = GL.glCreateProgram()
        GL.glAttachShader(program, my_vertex_shader)
        GL.glAttachShader(program, my_fragment_shader)
        if geometry is not None:
            GL.glAttachShader(program, my_geometry_shader)
        GL.glLinkProgram(program)
        return program
        
    def _create_shader(self, shader_prog, shader_type):
        """ Creates, links to a source, compiles and returns a shader.
            
            Keyword arguments:
            shader -- The shader text to use
            shader_type -- The OpenGL enum type of shader, it can be:
                           GL.GL_VERTEX_SHADER, GL.GL_GEOMETRY_SHADER or
                           GL.GL_FRAGMENT_SHADER
            
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
    
    def _load_matrices(self):
        """ Function doc """
        model = GL.glGetUniformLocation(self.gl_program, 'model_mat')
        GL.glUniformMatrix4fv(model, 1, GL.GL_FALSE, self.model_mat)
        view = GL.glGetUniformLocation(self.gl_program, 'view_mat')
        GL.glUniformMatrix4fv(view, 1, GL.GL_FALSE, self.view_mat)
        proj = GL.glGetUniformLocation(self.gl_program, 'projection_mat')
        GL.glUniformMatrix4fv(proj, 1, GL.GL_FALSE, self.proj_mat)
    
    def _load_lights(self):
        """ Function doc
        """
        light_pos = GL.glGetUniformLocation(self.gl_program, 'my_light.position')
        GL.glUniform3fv(light_pos, 1, self.light_position)
        #light_col = GL.glGetUniformLocation(self.gl_program, 'my_light.color')
        #GL.glUniform3fv(light_col, 1, self.light_color)
        amb_coef = GL.glGetUniformLocation(self.gl_program, 'my_light.ambient_coef')
        GL.glUniform1fv(amb_coef, 1, self.light_ambient_coef)
        shiny = GL.glGetUniformLocation(self.gl_program, 'my_light.shininess')
        GL.glUniform1fv(shiny, 1, self.light_shininess)
        intensity = GL.glGetUniformLocation(self.gl_program, 'my_light.intensity')
        GL.glUniform3fv(intensity, 1, self.light_intensity)
        #spec_col = GL.glGetUniformLocation(self.gl_program, 'my_light.specular_color')
        #GL.glUniform3fv(spec_col, 1, self.light_specular_color)
        return True
    
    def render(self, area, context):
        """ Function doc """
        if not self.gl_programs_compiled:
            self._create_gl_programs()
            self.gl_programs_compiled = True
        GL.glClearColor(self.bckgrnd_color[0], self.bckgrnd_color[1], self.bckgrnd_color[2], self.bckgrnd_color[3])
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        
        if self.draw_type is None:
            print("No draw type specified, closing program.")
            quit()
        elif self.vao is None:
            if self.gl_data is None:
                print("No data present, please load data")
            else:
                self.make_vaos()
        elif self.draw_type == "circles":
            # GL.glEnable(GL.GL_DEPTH_TEST)
            GL.glUseProgram(self.gl_program)
            GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
            self._load_matrices()
            GL.glBindVertexArray(self.vao)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, self.gl_data.xyz.itemsize*self.gl_data.xyz.flatten().shape[0], self.gl_data.xyz.flatten(), GL.GL_STREAM_DRAW)
            GL.glDrawArrays(GL.GL_POINTS, 0, self.elements)
            GL.glDisable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
            # GL.glDisable(GL.GL_DEPTH_TEST)
            GL.glBindVertexArray(0)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
            GL.glUseProgram(0)
    
    def make_vaos(self):
        """ Function doc """
        self.vao, self.elements, self.coord_vbo = vaos.circles(self.gl_program, self.gl_data)
    
    def load_data(self, data):
        """ Function doc """
        self.gl_data = data
    
    def get_viewport_pos(self, x, y):
        """ Function doc """
        px = (2.0*x - self.width) / self.width
        py = (self.height - 2.0*y) / self.height
        return [px, py, self.z_near]
    
    def mouse_pressed(self, widget, event):
        """ Function doc """
        return True
        left = event.button == 1
        middle = event.button == 2
        right = event.button == 3
        self.mouse_rotate = left and not (middle or right)
        self.mouse_zoom = right and not (middle or left)
        self.mouse_pan = middle and not (right  or left)
        self.mouse_x = float(event.x)
        self.mouse_y = float(event.y)
        self.drag_pos_x, self.drag_pos_y, self.drag_pos_z = self.get_viewport_pos(self.mouse_x, self.mouse_y)
        self.dragging = False
        if event.button == 1:
            self.dx = 0.0
            self.dy = 0.0
        if event.button == 2:
            pass
        if event.button == 3:
            pass
        return True
    
    def mouse_released(self, widget, event):
        """ Function doc """
        return True
        self.mouse_rotate = False
        self.mouse_zoom = False
        self.mouse_pan = False
        if event.button == 1:
            pass
        if event.button == 2:
            pass
        if event.button == 3:
            pass
        return True
    
    def mouse_motion(self, widget, event):
        """ Function doc """
        return True
        x = float(event.x)
        y = float(event.y)
        dx = x - self.mouse_x
        dy = y - self.mouse_y
        if (dx==0 and dy==0):
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
            self._update_cam_pos()
            self.dx = dx
            self.dy = dy
            self.dragging = True
            self.queue_draw()
        return True
    
    def mouse_scroll(self, widget, event):
        """ Function doc """
        return True
        if self.ctrl:
            if event.direction == Gdk.ScrollDirection.UP:
                self.model_mat = cam.my_glTranslatef(self.model_mat, [0.0, 0.0, -self.scroll])
            if event.direction == Gdk.ScrollDirection.DOWN:
                self.model_mat = cam.my_glTranslatef(self.model_mat, [0.0, 0.0, self.scroll])
        else:
            pos_z = self.cam_pos[2]
            if event.direction == Gdk.ScrollDirection.UP:
                self.z_near -= self.scroll
                self.z_far += self.scroll
            if event.direction == Gdk.ScrollDirection.DOWN:
                if (self.z_far-self.scroll) >= (self.min_zfar):
                    if (self.z_far-self.scroll) > (self.z_near+self.scroll):
                        self.z_near += self.scroll
                        self.z_far -= self.scroll
            if (self.z_near >= self.min_znear):
                self.proj_mat = cam.my_glPerspectivef(self.fov, self.var, self.z_near, self.z_far)
            else:
                if self.z_far < (self.min_zfar+self.min_znear):
                    self.z_near -= self.scroll
                    self.z_far = self.min_znear + self.min_zfar
                self.proj_mat = cam.my_glPerspectivef(self.fov, self.var, self.min_znear, self.z_far)
        self.queue_draw()
    
    def _rotate_view(self, dx, dy, x, y):
        """ Function doc """
        return True
        angle = math.sqrt(dx**2+dy**2)/float(self.width+1)*36.0*np.linalg.norm(self.cam_pos)
        if self.ctrl:
            if abs(dx) >= abs(dy):
                if (y-self.height/2.0) < 0:
                    rot_mat = cam.my_glRotatef(np.identity(4), angle, [0.0, 0.0, dx])
                else:
                    rot_mat = cam.my_glRotatef(np.identity(4), angle, [0.0, 0.0, -dx])
            else:
                if (x-self.width/2.0) < 0:
                    rot_mat = cam.my_glRotatef(np.identity(4), angle, [0.0, 0.0, -dy])
                else:
                    rot_mat = cam.my_glRotatef(np.identity(4), angle, [0.0, 0.0, dy])
        else:
            rot_mat = cam.my_glRotatef(np.identity(4), angle, [-dy, -dx, 0.0])
        self.model_mat = cam.my_glMultiplyMatricesf(self.model_mat, rot_mat)
        return True
    
    def _pan_view(self, x, y):
        """ Function doc """
        px, py, pz = self.get_viewport_pos(x, y)
        pan_mat = cam.my_glTranslatef(np.identity(4, dtype=float),
            [(px-self.drag_pos_x)*self.z_far/10.0, 
             (py-self.drag_pos_y)*self.z_far/10.0, 
             (pz-self.drag_pos_z)*self.z_far/10.0])
        self.model_mat = cam.my_glMultiplyMatricesf(self.model_mat, pan_mat)
        self.drag_pos_x = px
        self.drag_pos_y = py
        self.drag_pos_z = pz
        return True
    
    def _zoom_view(self, dy):
        """ Function doc """
        return True
        delta = (((self.z_far-self.z_near)/2.0)+self.z_near)/200.0
        move_z = dy * delta
        moved_mat = cam.my_glTranslatef(self.view_mat, [0.0, 0.0, move_z])
        moved_pos = cam.get_xyz_coords(moved_mat)
        if moved_pos[2] > 0.101:
            self.view_mat = moved_mat
            self.z_near -= move_z
            self.z_far -= move_z
            if self.z_near >= self.min_znear:
                self.proj_mat = cam.my_glPerspectivef(self.fov, self.var, self.z_near, self.z_far)
            else:
                if self.z_far < (self.min_zfar+self.min_znear):
                    self.z_near += move_z
                    self.z_far = self.min_zfar+self.min_znear
                self.proj_mat = cam.my_glPerspectivef(self.fov, self.var, self.min_znear, self.z_far)
        else:
            pass
        return True
    
    def key_pressed(self, widget, event):
        """ Function doc """
        k_name = Gdk.keyval_name(event.keyval)
        func = getattr(self, '_pressed_' + k_name, None)
        if func:
            func()
        return True
    
    def key_released(self, widget, event):
        """ Function doc """
        k_name = Gdk.keyval_name(event.keyval)
        func = getattr(self, '_released_' + k_name, None)
        if func:
            func()
        return True
    
    def _pressed_Escape(self):
        Gtk.main_quit()
    
    def _pressed_p(self):
        print(self.gl_data.xyz)
        print(self.gl_data.colors)
        print(self.gl_data.direc)
        # print(self.proj_mat)
        # print(self.model_mat)
        # print(self.view_mat)
        # print(self.get_cam_pos())
    
    def _pressed_r(self):
        import time
        self.animation = True
        while self.animation:
            self.gl_data.update_pos(1)
            self.get_window().invalidate_rect(None, False)
            self.get_window().process_updates(False)
            # time.sleep(.05)
        self.queue_render()
    
    def _pressed_s(self):
        import time
        self.animation = False
        self.queue_render()
    
