#!/usr/bin/env python3
# -*- coding: utf-8 -*-


vertex_shader_dot_simple = """
# version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_coord = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    frag_color = vert_color;
}
"""

fragment_shader_dot_simple = """
# version 330

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;

in vec3 frag_coord;
in vec3 frag_color;
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

sel_fragment_shader_dot_simple = """
# version 330

in vec3 frag_color;
out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""


vertex_shader_dot_circle = """
# version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_coord = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    frag_color = vert_color;
}
"""

fragment_shader_dot_circle = """
# version 330

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;

in vec3 frag_coord;
in vec3 frag_color;
out vec4 final_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    
    dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = (fog_end-dist)/(fog_end-fog_start);
        final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
    }
    else{
       final_color = vec4(frag_color, 1.0);
    }
}
"""

sel_fragment_shader_dot_circle = """
# version 330

in vec3 frag_color;
out vec4 final_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    final_color = vec4(frag_color, 1.0);
}
"""


vertex_shader_dot_disc = """
# version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_coord = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    frag_color = vert_color;
}
"""

fragment_shader_dot_disc = """
# version 330

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;

in vec3 frag_coord;
in vec3 frag_color;
out vec4 final_color;

void main(){
    vec4 circle_color;
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    float ligth_dist = length(gl_PointCoord - vec2(0.3, 0.3));
    circle_color = mix(vec4(frag_color, 1), vec4(0, 0, 0, 1), sqrt(ligth_dist)*.78);
    
    dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = (fog_end-dist)/(fog_end-fog_start);
        final_color = mix(fog_color, circle_color, fog_factor);
    }
    else{
       final_color = circle_color;
    }
}
"""

sel_fragment_shader_dot_disc = """
# version 330

in vec3 frag_color;
out vec4 final_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    final_color = vec4(frag_color, 1.0);
}
"""


vertex_shader_dot_extra = """
# version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec4 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_color = vert_color;
}
"""

fragment_shader_dot_extra = """
# version 330
uniform mat4 proj_mat;

in vec4 frag_coord;
in vec3 frag_color;
out vec4 final_color;

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
    vec4 pos = frag_coord;
    pos.z += z;
    pos = proj_mat * pos;
    gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    float ligth_dist = length(gl_PointCoord - vec2(0.3, 0.3));
    final_color = mix(vec4(frag_color, 1), vec4(0, 0, 0, 1), sqrt(ligth_dist)*.78);
}
"""

sel_fragment_shader_dot_extra = """
# version 330

in vec3 frag_color;
out vec4 final_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    final_color = vec4(frag_color, 1.0);
}
"""


shader_type = {0: { "vertex_shader"      : vertex_shader_dot_simple,
                    "fragment_shader"    : fragment_shader_dot_simple,
                    "sel_vertex_shader"  : vertex_shader_dot_simple,
                    "sel_fragment_shader": sel_fragment_shader_dot_simple
                  },
               1: {"vertex_shader"      : vertex_shader_dot_circle,
                   "fragment_shader"    : fragment_shader_dot_circle ,
                   "sel_vertex_shader"  : vertex_shader_dot_circle,
                   "sel_fragment_shader": sel_fragment_shader_dot_circle
                   },
               2: {"vertex_shader"      : vertex_shader_dot_disc,
                   "fragment_shader"    : fragment_shader_dot_disc,
                   "sel_vertex_shader"  : vertex_shader_dot_disc,
                   "sel_fragment_shader": sel_fragment_shader_dot_disc
                   },
               3: {"vertex_shader"      : vertex_shader_dot_extra,
                   "fragment_shader"    : fragment_shader_dot_extra,
                   "sel_vertex_shader"  : vertex_shader_dot_extra,
                   "sel_fragment_shader": sel_fragment_shader_dot_extra
                   },
}
