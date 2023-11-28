#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#


from __future__ import annotations

__version__ = "0.0.1"
__author__ = "Carlos Eduardo Sequeiros-Borja, José Fernando Ruggiero Bachega"
__mail__ = "casebor@gmail.com, ferbachega@gmail.com"


vertex_shader_lines = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

geometry_shader_lines = """
#version 330
precision highp float; 
precision highp int;

layout (lines) in;
layout (line_strip, max_vertices = 4) out;

uniform mat4 proj_mat;

in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;
out vec4 frag_coord;

void main(){
    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
    gl_Position = proj_mat * geom_coord[0];
    frag_color = geom_color[0];
    frag_coord = geom_coord[0];
    EmitVertex();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[0];
    frag_coord = mid_coord;
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[1];
    frag_coord = mid_coord;
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_coord = geom_coord[1];
    frag_color = geom_color[1];
    EmitVertex();
    EndPrimitive();
}
"""

fragment_shader_lines = """
#version 330
precision highp float; 
precision highp int;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;

in vec3 frag_color;
in vec4 frag_coord;

out vec4 final_color;

void main(){
    float dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = (fog_end-dist)/(fog_end-fog_start);
        final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
    }
    else{
       final_color = vec4(frag_color, 1.0);
    }
}
"""

################################## SELECTION ###################################

sel_vertex_shader_lines = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

sel_geometry_shader_lines = """
#version 330
precision highp float; 
precision highp int;

layout (lines) in;
layout (line_strip, max_vertices = 4) out;

uniform mat4 proj_mat;

in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;

void main(){
    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
    gl_Position = proj_mat * geom_coord[0];
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[0];
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_color = geom_color[1];
    EmitVertex();
    EndPrimitive();
}
"""

sel_fragment_shader_lines = """
#version 330
precision highp float; 
precision highp int;

in vec3 frag_color;

out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""


#################################### DOUBLE ####################################

vertex_shader_double = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

geometry_shader_double = """
#version 330
precision highp float; 
precision highp int;

layout (lines) in;
layout (line_strip, max_vertices = 8) out;

uniform mat4 proj_mat;
uniform mat4 view_mat;
const float separation = 0.04;

in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;
out vec4 frag_coord;

vec4 get_vec(vec3 point_A, vec3 point_B, vec3 campos){
    vec3 tmp = (point_B - point_A);
    tmp = cross(campos, tmp);
    tmp = normalize(tmp);
    return vec4(tmp, 0.0);
}

void main(){
    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
    mat4 inv_view = inverse(view_mat);
    vec3 campos = vec3(inv_view[3][0], inv_view[3][1], inv_view[3][2]);
    vec4 vort = get_vec(geom_coord[0].xyz, geom_coord[1].xyz, campos);
    
    gl_Position = proj_mat * (geom_coord[0] - vort*separation);
    frag_color = geom_color[0];
    frag_coord = geom_coord[0] - vort*separation;
    EmitVertex();
    gl_Position = proj_mat * (mid_coord - vort*separation);
    frag_color = geom_color[0];
    frag_coord = mid_coord - vort*separation;
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * (mid_coord - vort*separation);
    frag_color = geom_color[1];
    frag_coord = mid_coord - vort*separation;
    EmitVertex();
    gl_Position = proj_mat * (geom_coord[1] - vort*separation);
    frag_color = geom_color[1];
    frag_coord = geom_coord[1] - vort*separation;
    EmitVertex();
    EndPrimitive();
    
    gl_Position = proj_mat * (geom_coord[0] + vort*separation);
    frag_color = geom_color[0];
    frag_coord = geom_coord[0] + vort*separation;
    EmitVertex();
    gl_Position = proj_mat * (mid_coord + vort*separation);
    frag_color = geom_color[0];
    frag_coord = mid_coord + vort*separation;
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * (mid_coord + vort*separation);
    frag_color = geom_color[1];
    frag_coord = mid_coord + vort*separation;
    EmitVertex();
    gl_Position = proj_mat * (geom_coord[1] + vort*separation);
    frag_color = geom_color[1];
    frag_coord = geom_coord[1] + vort*separation;
    EmitVertex();
    EndPrimitive();
}
"""

fragment_shader_double = """
#version 330
precision highp float; 
precision highp int;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;

in vec3 frag_color;
in vec4 frag_coord;

out vec4 final_color;

void main(){
    float dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = (fog_end-dist)/(fog_end-fog_start);
        final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
    }
    else{
       final_color = vec4(frag_color, 1.0);
    }
}
"""


################################## SELECTION ###################################

sel_vertex_shader_double = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

sel_geometry_shader_double = """
#version 330
precision highp float; 
precision highp int;

layout (lines) in;
layout (line_strip, max_vertices = 8) out;

uniform mat4 proj_mat;
uniform mat4 view_mat;
const float separation = 0.04;

in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;

vec4 get_vec(vec3 point_A, vec3 point_B, vec3 campos){
    vec3 tmp = (point_B - point_A);
    tmp = cross(campos, tmp);
    tmp = normalize(tmp);
    return vec4(tmp, 0.0);
}

void main(){
    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
    mat4 inv_view = inverse(view_mat);
    vec3 campos = vec3(inv_view[3][0], inv_view[3][1], inv_view[3][2]);
    vec4 vort = get_vec(geom_coord[0].xyz, geom_coord[1].xyz, campos);
    
    gl_Position = proj_mat * (geom_coord[0] - vort*separation);
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * (mid_coord - vort*separation);
    frag_color = geom_color[0];
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * (mid_coord - vort*separation);
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * (geom_coord[1] - vort*separation);
    frag_color = geom_color[1];
    EmitVertex();
    EndPrimitive();
    
    gl_Position = proj_mat * (geom_coord[0] + vort*separation);
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * (mid_coord + vort*separation);
    frag_color = geom_color[0];
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * (mid_coord + vort*separation);
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * (geom_coord[1] + vort*separation);
    frag_color = geom_color[1];
    EmitVertex();
    EndPrimitive();
}
"""

sel_fragment_shader_double = """
#version 330
precision highp float; 
precision highp int;

in vec3 frag_color;

out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""



#################################### TRIPLE ####################################

vertex_shader_triple = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

geometry_shader_triple = """
#version 330
precision highp float; 
precision highp int;

layout (lines) in;
layout (line_strip, max_vertices = 12) out;

uniform mat4 proj_mat;
uniform mat4 view_mat;
const float separation = 0.06;

in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;
out vec4 frag_coord;

vec4 get_vec(vec3 point_A, vec3 point_B, vec3 campos){
    vec3 tmp = (point_B - point_A);
    tmp = cross(campos, tmp);
    tmp = normalize(tmp);
    return vec4(tmp, 0.0);
}

void main(){
    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
    mat4 inv_view = inverse(view_mat);
    vec3 campos = vec3(inv_view[3][0], inv_view[3][1], inv_view[3][2]);
    vec4 vort = get_vec(geom_coord[0].xyz, geom_coord[1].xyz, campos);
    
    gl_Position = proj_mat * (geom_coord[0] - vort*separation);
    frag_color = geom_color[0];
    frag_coord = geom_coord[0] - vort*separation;
    EmitVertex();
    gl_Position = proj_mat * (mid_coord - vort*separation);
    frag_color = geom_color[0];
    frag_coord = mid_coord - vort*separation;
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * (mid_coord - vort*separation);
    frag_color = geom_color[1];
    frag_coord = mid_coord - vort*separation;
    EmitVertex();
    gl_Position = proj_mat * (geom_coord[1] - vort*separation);
    frag_color = geom_color[1];
    frag_coord = geom_coord[1] - vort*separation;
    EmitVertex();
    EndPrimitive();
    
    gl_Position = proj_mat * geom_coord[0];
    frag_color = geom_color[0];
    frag_coord = geom_coord[0];
    EmitVertex();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[0];
    frag_coord = mid_coord;
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[1];
    frag_coord = mid_coord;
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_color = geom_color[1];
    frag_coord = geom_coord[1];
    EmitVertex();
    EndPrimitive();
    
    gl_Position = proj_mat * (geom_coord[0] + vort*separation);
    frag_color = geom_color[0];
    frag_coord = geom_coord[0] + vort*separation;
    EmitVertex();
    gl_Position = proj_mat * (mid_coord + vort*separation);
    frag_color = geom_color[0];
    frag_coord = mid_coord + vort*separation;
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * (mid_coord + vort*separation);
    frag_color = geom_color[1];
    frag_coord = mid_coord + vort*separation;
    EmitVertex();
    gl_Position = proj_mat * (geom_coord[1] + vort*separation);
    frag_color = geom_color[1];
    frag_coord = geom_coord[1] + vort*separation;
    EmitVertex();
    EndPrimitive();
}
"""

fragment_shader_triple = """
#version 330
precision highp float; 
precision highp int;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;

in vec3 frag_color;
in vec4 frag_coord;

out vec4 final_color;

void main(){
    float dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = (fog_end-dist)/(fog_end-fog_start);
        final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
    }
    else{
       final_color = vec4(frag_color, 1.0);
    }
}
"""


################################## SELECTION ###################################

sel_vertex_shader_triple = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

sel_geometry_shader_triple = """
#version 330
precision highp float; 
precision highp int;

layout (lines) in;
layout (line_strip, max_vertices = 12) out;

uniform mat4 proj_mat;
uniform mat4 view_mat;
const float separation = 0.06;

in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;

vec4 get_vec(vec3 point_A, vec3 point_B, vec3 campos){
    vec3 tmp = (point_B - point_A);
    tmp = cross(campos, tmp);
    tmp = normalize(tmp);
    return vec4(tmp, 0.0);
}

void main(){
    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
    mat4 inv_view = inverse(view_mat);
    vec3 campos = vec3(inv_view[3][0], inv_view[3][1], inv_view[3][2]);
    vec4 vort = get_vec(geom_coord[0].xyz, geom_coord[1].xyz, campos);
    
    gl_Position = proj_mat * (geom_coord[0] - vort*separation);
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * (mid_coord - vort*separation);
    frag_color = geom_color[0];
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * (mid_coord - vort*separation);
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * (geom_coord[1] - vort*separation);
    frag_color = geom_color[1];
    EmitVertex();
    EndPrimitive();
    
    gl_Position = proj_mat * geom_coord[0];
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[0];
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_color = geom_color[1];
    EmitVertex();
    EndPrimitive();
    
    gl_Position = proj_mat * (geom_coord[0] + vort*separation);
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * (mid_coord + vort*separation);
    frag_color = geom_color[0];
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * (mid_coord + vort*separation);
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * (geom_coord[1] + vort*separation);
    frag_color = geom_color[1];
    EmitVertex();
    EndPrimitive();
}
"""

sel_fragment_shader_triple = """
#version 330
precision highp float; 
precision highp int;

in vec3 frag_color;

out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""


shader_type ={
            0: { "vertex_shader"      : vertex_shader_lines,
                 "geometry_shader"    : geometry_shader_lines,
                 "fragment_shader"    : fragment_shader_lines,
                 "sel_vertex_shader"  : sel_vertex_shader_lines,
                 "sel_geometry_shader": sel_geometry_shader_lines,
                 "sel_fragment_shader": sel_fragment_shader_lines
               },
            1: {"vertex_shader"      : vertex_shader_double,
                "geometry_shader"    : geometry_shader_double,
                "fragment_shader"    : fragment_shader_double,
                "sel_vertex_shader"  : sel_vertex_shader_double,
                "sel_geometry_shader": sel_geometry_shader_double,
                "sel_fragment_shader": sel_fragment_shader_double
                },
            2: {"vertex_shader"      : vertex_shader_triple,
                "geometry_shader"    : geometry_shader_triple,
                "fragment_shader"    : fragment_shader_triple,
                "sel_vertex_shader"  : sel_vertex_shader_triple,
                "sel_geometry_shader": sel_geometry_shader_triple,
                "sel_fragment_shader": sel_fragment_shader_triple
                }
}
