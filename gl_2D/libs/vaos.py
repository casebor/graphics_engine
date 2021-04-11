#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import ctypes
import numpy as np

from OpenGL import GL

def circles(program, data):
    """ Function doc """
    vertex_array_object = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vertex_array_object)
    coords = data.xyz.flatten()
    colors = data.colors.flatten()
    radii = data.radii
    
    coord_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*coords.shape[0], coords, GL.GL_STREAM_DRAW)
    gl_coords = GL.glGetAttribLocation(program, 'vert_coord')
    GL.glEnableVertexAttribArray(gl_coords)
    GL.glVertexAttribPointer(gl_coords, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
    
    col_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.itemsize*colors.shape[0], colors, GL.GL_STATIC_DRAW)
    gl_colors = GL.glGetAttribLocation(program, 'vert_color')
    GL.glEnableVertexAttribArray(gl_colors)
    GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
    
    radii_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, radii_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, radii.itemsize*radii.shape[0], radii, GL.GL_STATIC_DRAW)
    gl_radii = GL.glGetAttribLocation(program, 'vert_radius')
    GL.glEnableVertexAttribArray(gl_radii)
    GL.glVertexAttribPointer(gl_radii, 1, GL.GL_FLOAT, GL.GL_FALSE, 1*radii.itemsize, ctypes.c_void_p(0))
    
    GL.glBindVertexArray(0)
    GL.glDisableVertexAttribArray(gl_coords)
    GL.glDisableVertexAttribArray(gl_colors)
    GL.glDisableVertexAttribArray(gl_radii)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
    return vertex_array_object, data.xyz.shape[0], coord_vbo

def _make_vbos(program, coords, colors, indexes=None, radii=None):
    """ Function doc """
    if indexes is not None:
        ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.itemsize*indexes.shape[0], indexes, GL.GL_DYNAMIC_DRAW)
    
    coord_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*coords.shape[0], coords, GL.GL_STATIC_DRAW)
    gl_coords = GL.glGetAttribLocation(program, 'vert_coord')
    GL.glEnableVertexAttribArray(gl_coords)
    GL.glVertexAttribPointer(gl_coords, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
    
    col_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.itemsize*colors.shape[0], colors, GL.GL_STATIC_DRAW)
    gl_colors = GL.glGetAttribLocation(program, 'vert_color')
    GL.glEnableVertexAttribArray(gl_colors)
    GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
    
    if radii is not None:
        radii_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, radii_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, radii.itemsize*radii.shape[0], radii, GL.GL_STATIC_DRAW)
        gl_radii = GL.glGetAttribLocation(program, 'vert_radius')
        GL.glEnableVertexAttribArray(gl_radii)
        GL.glVertexAttribPointer(gl_radii, 1, GL.GL_FLOAT, GL.GL_FALSE, 1*radii.itemsize, ctypes.c_void_p(0))
        return gl_coords, gl_colors, gl_radii
    
    return gl_coords, gl_colors

def points(program, data):
    """ Function doc """
    vertex_array_object = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vertex_array_object)
    coords = data.xyz.flatten()
    colors = data.colors.flatten()
    rads = None
    if data.radii is not None:
        pos, cols, rads = _make_vbos(program, coords, colors, radii=data.radii)
    else:
        pos, cols = _make_vbos(program, coords, colors)
    
    GL.glBindVertexArray(0)
    GL.glDisableVertexAttribArray(pos)
    GL.glDisableVertexAttribArray(cols)
    if rads is not None:
        GL.glDisableVertexAttribArray(rads)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
    return vertex_array_object, data.xyz.shape[0]

def lines(program, data):
    """ Function doc """
    vertex_array_object = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vertex_array_object)
    coords = data.xyz.flatten()
    colors = data.colors.flatten()
    indexes = data.indexes
    pos, cols = _make_vbos(program, coords, colors, indexes)
    
    GL.glBindVertexArray(0)
    GL.glDisableVertexAttribArray(pos)
    GL.glDisableVertexAttribArray(cols)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
    return vertex_array_object, indexes.shape[0]

def triangles(program, data):
    vertex_array_object = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vertex_array_object)
    coords = data.xyz.flatten()
    colors = data.colors.flatten()
    normals = data.normals.flatten()
    indexes = data.indexes
    pos, cols, norms = _make_vbos(program, coords, colors, indexes, normals)
    
    GL.glBindVertexArray(0)
    GL.glDisableVertexAttribArray(pos)
    GL.glDisableVertexAttribArray(cols)
    GL.glDisableVertexAttribArray(norms)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
    return vertex_array_object, indexes.shape[0]




def make_text_texture():
    """ Function doc """
    from PIL import Image
    image_a = Image.open("SimSun_ExtB.tga")
    #print('opened file: size=', image_a.size, 'format=', image_a.format)
    ix = image_a.size[0]
    iy = image_a.size[1]
    image_a = np.array(list(image_a.getdata()), np.uint8)
    text_texture = GL.glGenTextures(1)
    GL.glActiveTexture(GL.GL_TEXTURE0)
    GL.glBindTexture(GL.GL_TEXTURE_2D, text_texture)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_REPEAT)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_REPEAT)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
    GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, ix, iy, 0, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, image_a)
    return text_texture

def make_text(program):
    """ Function doc """
    phrase = "Hello World!!!"
    text_id = np.zeros(len(phrase),np.uint32)
    indexes = np.zeros(len(phrase),np.uint32)
    for i,letter in enumerate(phrase):
        text_id[i] = ord(letter)
        indexes[i] = i
    coords = np.zeros(len(phrase)*3,np.float32)
    point = [-1, 1, 0]
    for i in range(0, coords.size, 3):
        coords[i] = point[0] + i*0.06
        coords[i+1] = point[1]
        coords[i+2] = point[2]
    
    #coords = np.array([-1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
    #                   -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    #                   -1.0,-1.0, 0.0, 0.0,-1.0, 0.0, 1.0,-1.0, 0.0],dtype=np.float32)
    #text_id = np.array([33, 34, 35, 56, 111, 122, 87, 90, 666], dtype=np.int32)
    #indexes = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=np.uint32)
    
    vertex_array_object = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vertex_array_object)
    
    ind_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
    GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.itemsize*int(len(indexes)), indexes, GL.GL_DYNAMIC_DRAW)
    
    coord_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
    gl_coord = GL.glGetAttribLocation(program, 'vert_coord')
    GL.glEnableVertexAttribArray(gl_coord)
    GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
    
    tex_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, tex_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, text_id.itemsize*len(text_id), text_id, GL.GL_STATIC_DRAW)
    gl_texture = GL.glGetAttribLocation(program, 'vert_id')
    GL.glEnableVertexAttribArray(gl_texture)
    GL.glVertexAttribPointer(gl_texture, 1, GL.GL_FLOAT, GL.GL_FALSE, 1*text_id.itemsize, ctypes.c_void_p(0))
    
    GL.glBindVertexArray(0)
    GL.glDisableVertexAttribArray(gl_coord)
    GL.glDisableVertexAttribArray(gl_texture)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
    return vertex_array_object, (ind_vbo, coord_vbo, tex_vbo), int(len(indexes))

def make_texture_texture():
    """ Function doc """
    from PIL import Image
    image_a = Image.open("test.tga")
    #print('opened file: size=', image_a.size, 'format=', image_a.format)
    ix = image_a.size[0]
    iy = image_a.size[1]
    image_a = np.array(list(image_a.getdata()),np.uint8)
    tex_texture = GL.glGenTextures(1)
    #tex_texture = GL.glGenTextures(2)
    GL.glActiveTexture(GL.GL_TEXTURE0)
    GL.glBindTexture(GL.GL_TEXTURE_2D, tex_texture)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_REPEAT)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_REPEAT)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
    GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, ix, iy, 0, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, image_a)
    #image_b = img.open("cry.bmp")
    #ix = image_b.size[0]
    #iy = image_b.size[1]
    #image_b = image_b.tobytes("raw", "RGBX", 0, -1)
    #GL.glActiveTexture(GL.GL_TEXTURE1)
    #GL.glBindTexture(GL.GL_TEXTURE_2D, tex_texture[1])
    #GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_REPEAT)
    #GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_REPEAT)
    #GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
    #GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
    #GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, ix, iy, 0, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, image_b)
    return tex_texture

def make_texture(program):
    """ Function doc """
    coords = np.array([-1.0, 1.0, 0.0,-1.0,-1.0, 0.0, 1.0,-1.0, 0.0,
                       -1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0,-1.0, 0.0,],dtype=np.float32)
    textur = np.array([ 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                        0.0, 1.0, 1.0, 1.0, 1.0, 0.0],dtype=np.float32)
    
    vertex_array_object = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vertex_array_object)
    
    coord_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
    gl_coord = GL.glGetAttribLocation(program, 'vert_coord')
    GL.glEnableVertexAttribArray(gl_coord)
    GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
    
    tex_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, tex_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, textur.itemsize*len(textur), textur, GL.GL_STATIC_DRAW)
    gl_texture = GL.glGetAttribLocation(program, 'vert_text')
    GL.glEnableVertexAttribArray(gl_texture)
    GL.glVertexAttribPointer(gl_texture, 2, GL.GL_FLOAT, GL.GL_FALSE, 2*textur.itemsize, ctypes.c_void_p(0))
    
    GL.glBindVertexArray(0)
    GL.glDisableVertexAttribArray(gl_coord)
    GL.glDisableVertexAttribArray(gl_texture)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
    return vertex_array_object, (coord_vbo, tex_vbo), int(len(coords)/3)

def make_dots_surface(program):
    """ Function doc """
    vertex_array_object = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vertex_array_object)
    coords = np.array([-1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                       -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       -1.0,-1.0, 0.0, 0.0,-1.0, 0.0, 1.0,-1.0, 0.0],dtype=np.float32)
    colors = np.array([ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
                        1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0,
                        0.5, 0.5, 0.5, 0.2, 0.3, 0.4, 0.9, 0.5, 0.1],dtype=np.float32)
    indexes = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=np.uint32)
    
    ind_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
    GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.itemsize*int(len(indexes)), indexes, GL.GL_DYNAMIC_DRAW)
    
    coord_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
    gl_coord = GL.glGetAttribLocation(program, 'vert_coord')
    GL.glEnableVertexAttribArray(gl_coord)
    GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
    
    col_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.itemsize*len(colors), colors, GL.GL_STATIC_DRAW)
    gl_color = GL.glGetAttribLocation(program, 'vert_color')
    GL.glEnableVertexAttribArray(gl_color)
    GL.glVertexAttribPointer(gl_color, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
    
    GL.glBindVertexArray(0)
    GL.glDisableVertexAttribArray(gl_coord)
    GL.glDisableVertexAttribArray(gl_color)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
    return vertex_array_object, (ind_vbo, coord_vbo, col_vbo), int(len(indexes))

def make_icosahedron(program):
    """ Function doc """
    vertex_array_object = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vertex_array_object)
    coords = np.array([-1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                       -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       -1.0,-1.0, 0.0, 0.0,-1.0, 0.0, 1.0,-1.0, 0.0],dtype=np.float32)
    colors = np.array([ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
                        1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0,
                        0.5, 0.5, 0.5, 0.2, 0.3, 0.4, 0.9, 0.5, 0.1],dtype=np.float32)
    indexes = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=np.uint32)
    
    ind_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
    GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.itemsize*int(len(indexes)), indexes, GL.GL_DYNAMIC_DRAW)
    
    coord_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
    gl_coord = GL.glGetAttribLocation(program, 'vert_coord')
    GL.glEnableVertexAttribArray(gl_coord)
    GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
    
    col_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.itemsize*len(colors), colors, GL.GL_STATIC_DRAW)
    gl_color = GL.glGetAttribLocation(program, 'vert_color')
    GL.glEnableVertexAttribArray(gl_color)
    GL.glVertexAttribPointer(gl_color, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
    
    GL.glBindVertexArray(0)
    GL.glDisableVertexAttribArray(gl_coord)
    GL.glDisableVertexAttribArray(gl_color)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
    return vertex_array_object, (ind_vbo, coord_vbo, col_vbo), int(len(indexes))

def make_sphere(program, level='level_2'):
    """ Function doc """
    nucleus = [-1.0, 1.0, 0.0, 1.0, 1.0, 0.0,-1.0,-1.0, 0.0, 1.0,-1.0, 0.0]
    colores = [ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5]
    qtty = int(len(nucleus)/3)
    coords = np.array([], dtype=np.float32)
    centers = np.array([], dtype=np.float32)
    colors = np.array([], dtype=np.float32)
    indexes = np.array([], dtype=np.uint32)
    for i in range(qtty):
        crds = np.copy(sphd.sphere_vertices[level])
        inds = np.copy(sphd.sphere_triangles[level])
        offset = int(len(crds)/3)
        cols = np.array(colores[i*3:(i+1)*3]*offset, dtype=np.float32)
        cnts = np.array(nucleus[i*3:(i+1)*3]*offset, dtype=np.float32)
        for j in range(offset):
            crds[j*3] = crds[j*3] + nucleus[i*3]
            crds[j*3+1] = crds[j*3+1] + nucleus[i*3+1]
            crds[j*3+2] = crds[j*3+2] + nucleus[i*3+2]
        inds += i*offset
        coords = np.concatenate((coords, crds))
        centers = np.concatenate((centers, cnts))
        colors = np.concatenate((colors, cols))
        indexes = np.concatenate((indexes, inds))
    #coords = np.array(coords, dtype=np.float32)
    #indexes = np.array(indexes, dtype=np.uint32)
    #colors = np.array(colors, dtype=np.float32)
    
    vertex_array_object = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vertex_array_object)
    
    ind_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
    GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.itemsize*int(len(indexes)), indexes, GL.GL_DYNAMIC_DRAW)
    
    coord_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
    gl_coord = GL.glGetAttribLocation(program, 'vert_coord')
    GL.glEnableVertexAttribArray(gl_coord)
    GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
    
    centr_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, centr_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, centers.itemsize*len(centers), centers, GL.GL_STATIC_DRAW)
    gl_center = GL.glGetAttribLocation(program, 'vert_centr')
    GL.glEnableVertexAttribArray(gl_center)
    GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*centers.itemsize, ctypes.c_void_p(0))
    
    col_vbo = GL.glGenBuffers(1)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
    GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.itemsize*len(colors), colors, GL.GL_STATIC_DRAW)
    gl_colors = GL.glGetAttribLocation(program, 'vert_color')
    GL.glEnableVertexAttribArray(gl_colors)
    GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
    
    GL.glBindVertexArray(0)
    GL.glDisableVertexAttribArray(gl_coord)
    GL.glDisableVertexAttribArray(gl_center)
    GL.glDisableVertexAttribArray(gl_colors)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
    return vertex_array_object, (ind_vbo, coord_vbo, col_vbo), int(len(indexes))
