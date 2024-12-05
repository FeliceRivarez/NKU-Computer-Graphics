#pragma once
#ifndef __CONTAINER_H__
#define __CONTAINER_H__

#include <glad/glad.h>
#include <glfw3.h>
#include "Shader.h"
#include "Configure.h"
#include "Camera.h"

namespace Glb {
	class Container {
	public:
		void resetSize(float x, float y, float z);
		void draw();
		void init();

	private:

		Shader* shader = NULL;

		float x = 1.0f;
		float y = 1.0f;
		float z = 1.0f;

		GLuint VBO = 0;
		GLuint VAO = 0;

	};
}
#endif