#version 330 core
out vec4 FragColor;
in vec2 texCoord;

uniform sampler2D aTexture;
uniform float contrast;

void main()
{
	FragColor = texture(aTexture, texCoord) * contrast;
}