#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  dots.py
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

vertex_shader_dot_simple = """
# version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 v_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    v_color     = vert_color; 
}
"""

fragment_shader_dot_simple = """
# version 330

uniform vec4 fog_color;
uniform float fog_start;

in vec3 v_color;
out vec4 out_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    float dist2 = length(gl_PointCoord.xy - vec2(0.4,0.4));
    float sphere_factor = pow( 1.025 - dist2 , 1.75 ); 
    out_color = vec4(v_color * sphere_factor, 1.0);
}
"""

vertex_shader_dot_simple2 = """
# version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 v_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    v_color     = vert_color; 
}
"""

fragment_shader_dot_simple2 = """
# version 330

in vec3 v_color;
out vec4 out_color;

void main()
{
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    float ligth_dist = length(gl_PointCoord - vec2(0.3, 0.3));
    out_color = mix(vec4(v_color, 1), vec4(0, 0, 0, 1), sqrt(ligth_dist)*.78);
}
"""

vertex_shader_dot_simple3 = """
# version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 v_color;
out vec4 frag_pos;

void main(){
    frag_pos = view_mat * model_mat * vec4(vert_coord, 1.0);
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    v_color     = vert_color; 
}
"""

fragment_shader_dot_simple3 = """
# version 330
uniform mat4 proj_mat;

in vec3 v_color;
in vec4 frag_pos;
out vec4 out_color;

void main(){
    vec2 P = gl_PointCoord.xy - vec2(0.5,0.5);
    float point_size = 5.0;
    float distance = length(P);
    vec2 texcoord = gl_PointCoord* 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;
    if (d <= 0.0)
        discard;
    float z = sqrt(d);
    vec4 pos = frag_pos;
    pos.z += z;
    pos = proj_mat * pos;
    gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    float ligth_dist = length(gl_PointCoord - vec2(0.3, 0.3));
    out_color = mix(vec4(v_color, 1), vec4(0, 0, 0, 1), sqrt(ligth_dist)*.78);
}
"""


shader_type = {0: { "vertex_shader"      : vertex_shader_dot_simple,
                    "fragment_shader"    : fragment_shader_dot_simple,
                    "sel_vertex_shader"  : vertex_shader_dot_simple,
                    "sel_fragment_shader": fragment_shader_dot_simple
                  },
               1: {"vertex_shader"      : vertex_shader_dot_simple2,
                   "fragment_shader"    : fragment_shader_dot_simple2 ,
                   "sel_vertex_shader"  : vertex_shader_dot_simple,
                   "sel_fragment_shader": fragment_shader_dot_simple
                   },
               2: {"vertex_shader"      : vertex_shader_dot_simple3,
                   "fragment_shader"    : fragment_shader_dot_simple3,
                   "sel_vertex_shader"  : vertex_shader_dot_simple,
                   "sel_fragment_shader": fragment_shader_dot_simple
                   }
}
