
#ifdef GL_ES
precision mediump float;
#endif

uniform vec2    u_resolution;
uniform vec2    u_mouse;
uniform float   u_time;
varying vec4    v_position;
varying vec3    v_normal;
varying vec4    v_color;

float sphIntersect(vec3 ro, vec3 rd, vec3 sph, float rad){
    vec3 oc = ro - sph;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - rad*rad;
    float h = b*b - c;
    if( h<0.0 ) return -1.0;
    return -b - sqrt(h);
}

void main(void) {
    float right = u_resolution.x / u_resolution.y;
    float left = -right;
    vec2 uBottomLeft = vec2(left, -1.0);
    vec2 uTopRight = vec2(right, 1.0);
    
    vec3 color = vec3(1.0);
    vec2 st = gl_FragCoord.xy/u_resolution.xy;
    st.x *= u_resolution.x/u_resolution.y;
    color = vec3(st.x,st.y,abs(sin(u_time)));
    
    vec3 r0 = vec3(uBottomLeft + (gl_FragCoord.xy/u_resolution) * (uTopRight - uBottomLeft), 1.0);
    //vec3 r0 = vec3(gl_FragCoord.xy/u_resolution, 1.0);
    vec3 rd = vec3(0.0, 0.0, -1.0);
    //float t = raySphereIntersect(r0, rd);
    float t = sphIntersect(r0, rd, vec3(v_normal.xy, 0.0), 0.6);
    if (t < 0.0) {
        discard;
    }
    vec3 coord = r0 + rd * t;
    vec3 normal = normalize(coord - v_normal.xyz);
    gl_FragColor = vec4(v_color.xyz, 1.0);
    //gl_FragColor = vec4(color, 1.0);
}
