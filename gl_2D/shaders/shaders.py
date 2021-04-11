#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

vertex_shader = """
#version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 projection_mat;
const vec3 ligth_pos = vec3(1,1,2);

in vec3 vert_coord;
in vec3 vert_color;
in float vert_radius;

out vec3 frag_color;
//out vec4 frag_coord;
out vec3 v_light_direction;

void main(){
    frag_color = vert_color;
    vec4 frag_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
    v_light_direction = normalize(ligth_pos);
    gl_Position = projection_mat * frag_coord;
    gl_PointSize = vert_radius;
}
"""

fragment_shader = """
#version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 projection_mat;

in vec3 frag_color;
//in vec4 frag_coord;
in vec3 v_light_direction;

out vec4 final_color;

void main(){
    float dist = length(gl_PointCoord - vec2(0.5, 0.5));
    if (dist > 0.5)
        discard;
    float ligth_dist = length(gl_PointCoord - v_light_direction.xy);
    final_color = mix(vec4(frag_color, 1), vec4(0, 0, 0, 1), sqrt(ligth_dist)*.78);
}
"""

vertex_shader_glumpy = """
#version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 projection_mat;
const vec3 ligth_pos = vec3(0,1,2);

in vec3 vert_coord;     // attribute vec3 position;
in vec3 vert_color;     // attribute vec3 color;
in float vert_radius;   // attribute float radius;
//const float vert_dot_size = 50.5;  // attribute float radius;

out vec3 frag_color;    // varying vec3 v_color;
out float frag_radius;  // varying float v_radius;
out float frag_size;    // varying float v_size;
out vec4 frag_coord;    // varying vec4 v_eye_position;

varying vec3 v_light_direction;

void main (void){
    frag_color = vert_color;
    frag_radius = vert_radius;
    frag_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
    v_light_direction = normalize(ligth_pos);
    gl_Position = projection_mat * frag_coord;
    vec4 p = projection_mat * vec4(vert_radius, vert_radius, frag_coord.z, frag_coord.w);
    frag_size = 512.0 * p.x / p.w;
    gl_PointSize = frag_size + 5.0;
}
"""

fragment_shader_glumpy = """
#version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 projection_mat;

vec4 outline(float distance, float linewidth, float antialias, vec4 fg_color, vec4 bg_color){
    vec4 frag_color;
    float t = linewidth/2.0 - antialias;
    float signed_distance = distance;
    float border_distance = abs(signed_distance) - t;
    float alpha = border_distance/antialias;
    alpha = exp(-alpha*alpha);

    if( border_distance < 0.0 )
        frag_color = fg_color;
    else if( signed_distance < 0.0 )
        frag_color = mix(bg_color, fg_color, sqrt(alpha));
    else {
        if( abs(signed_distance) < (linewidth/2.0 + antialias) ) {
            frag_color = vec4(fg_color.rgb, fg_color.a * alpha);
        } else {
            discard;
        }
    }
    return frag_color;
}

in vec3 frag_color;       // varying vec3 v_color;
in float frag_radius;     // varying float v_radius;
in float frag_size;       // varying float v_size;
in vec4 frag_coord;       // varying vec4 v_eye_position;

varying vec3 v_light_direction;

void main(){
    vec2 P = gl_PointCoord.xy - vec2(0.5,0.5);
    float point_size = frag_size  + 5.0;
    float distance = length(P*point_size) - frag_size/2;
    vec2 texcoord = gl_PointCoord* 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;
    if (d <= 0.0) discard;
    float z = sqrt(d);
    vec4 pos = frag_coord;
    pos.z += frag_radius*z;
    vec3 pos2 = pos.xyz;
    pos = projection_mat * pos;
    gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    vec3 normal = vec3(x,y,z);
    float diffuse = clamp(dot(normal, v_light_direction), 0.0, 1.0);
    vec4 color = vec4((0.5 + 0.5*diffuse)*frag_color, 1.0);
    gl_FragColor = outline(distance, 1.0, 1.0, vec4(0,0,0,1), color);
}
"""

vertex_shader_triangles = """
#version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 projection_mat;

in vec3 vert_coord;
in vec3 vert_color;
in vec3 vert_norm;

out vec3 frag_coord;
out vec3 frag_color;
out vec3 frag_norm;

void main(){
    frag_coord = (view_mat * model_mat * vec4(vert_coord, 1)).xyz;
    frag_color = vert_color;
    //frag_norm = mat3(transpose(inverse(model_mat))) * vert_norm;
    frag_norm = vert_norm;
    gl_Position = projection_mat * vec4(frag_coord, 1);
}
"""

fragment_shader_triangles = """
#version 330

struct Light {
    vec3 position;
    //vec3 color;
    vec3 intensity;
    //vec3 specular_color;
    float ambient_coef;
    float shininess;
};

uniform Light my_light;

in vec3 frag_coord;
in vec3 frag_color;
in vec3 frag_norm;

out vec4 final_color;

void main(){
    vec3 normal = normalize(frag_norm);
    vec3 vert_to_light = normalize(my_light.position);
    vec3 vert_to_cam = normalize(frag_coord);
    
    // Ambient Component
    vec3 ambient = my_light.ambient_coef * frag_color * my_light.intensity;
    
    // Diffuse component
    float diffuse_coef = max(0.0, dot(normal, vert_to_light));
    vec3 diffuse = diffuse_coef * frag_color * my_light.intensity;
    
    // Specular component
    float specular_coef = 0.0;
    if (diffuse_coef > 0.0)
        specular_coef = pow(max(0.0, dot(vert_to_cam, reflect(-vert_to_light, normal))), my_light.shininess);
    vec3 specular = specular_coef * my_light.intensity;
    specular = specular * (vec3(1) - diffuse);
    
    final_color = vec4(ambient + diffuse + specular, 1.0);
}
"""