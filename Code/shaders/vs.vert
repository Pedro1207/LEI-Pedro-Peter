vs=`
#version 300 es

precision highp float;

in vec4 position;
in vec2 texcoord;

out vec2 texCoordV;

void main() {
    gl_Position = position;
    texCoordV = texcoord;
}
`