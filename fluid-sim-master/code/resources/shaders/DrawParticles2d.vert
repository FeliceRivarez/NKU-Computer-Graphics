
#version 330

layout(location = 0) in vec2 pointPosition;
layout(location = 1) in float density;

uniform float scale;

out float densNorm;

void main() {
    gl_PointSize = 3;
    gl_Position = vec4(pointPosition.xy/scale, 0.0, 1.0);

    densNorm = density / 5000.0f;   // normalization
}
