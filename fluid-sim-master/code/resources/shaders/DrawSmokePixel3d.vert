#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec2 vertTexCoord;

out vec2 texCoord;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
	gl_Position = projection * view * model * vec4(position, 1.0);
	texCoord = vec2(vertTexCoord.x, 1-vertTexCoord.y);
}