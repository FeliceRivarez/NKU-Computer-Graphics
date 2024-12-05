#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoord;
layout (location = 2) in float aColor;

out float ourColor;
out vec2 TexCoord;

void main()
{
	gl_Position = vec4(aPos, 0.0, 1.0);
	ourColor = aColor;
	TexCoord = vec2(aTexCoord.x, 1 - aTexCoord.y);
}