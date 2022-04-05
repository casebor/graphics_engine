#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#


vertex_shader_spheres = """
#version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

in vec3 vert_coord;
in vec3 vert_color;
in vec3 vert_instance;
in float vert_radius;

vec3 vert_norm;

out vec3 frag_coord;
out vec3 frag_color;
out vec3 frag_norm;

void main(){
    mat4 modelview = view_mat * model_mat;
    vec3 offset_coord = vert_coord * vert_radius + vert_instance;
    gl_Position = proj_mat * modelview * vec4(offset_coord, 1.0);
    
    vert_norm = normalize(offset_coord - vert_instance);
    frag_coord = vec3(modelview * vec4(offset_coord, 1.0));
    frag_norm = mat3(transpose(inverse(model_mat))) * vert_norm;
    frag_color = vert_color;
}
"""
fragment_shader_spheres = """
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

vec4 calculate_color(vec3 fnrm, vec3 fcrd, vec3 fcol){
    vec3 normal = normalize(fnrm);
    vec3 vert_to_light = normalize(my_light.position);
    vec3 vert_to_cam = normalize(fcrd);
    // Ambient Component
    vec3 ambient = my_light.ambient_coef * fcol * my_light.intensity;
    // Diffuse component
    float diffuse_coef = max(0.0, dot(normal, vert_to_light));
    vec3 diffuse = diffuse_coef * fcol * my_light.intensity;
    // Specular component
    float specular_coef = 0.0;
    if (diffuse_coef > 0.0)
        specular_coef = pow(max(0.0, dot(vert_to_cam, reflect(vert_to_light, normal))), my_light.shininess);
    vec3 specular = specular_coef * my_light.intensity;
    specular = specular * (vec3(1) - diffuse);
    vec4 out_color = vec4(ambient + diffuse + specular, 1.0);
    return out_color;
}

void main(){
    final_color = calculate_color(frag_norm, frag_coord, frag_color);
}
"""

