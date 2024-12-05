#pragma once
#ifndef __COMPONENT_H__
#define __COMPONENT_H__

#include <glfw3.h>

// 用于组织粒子/网格系统、渲染器、求解器

namespace Glb {
	class Component {
	public:
		int id;
		char* description;

		Component() {
			id = -1;
			this->description = NULL;
		}
		Component(char* description) {
			id = -1;
			this->description = description;
		}

		virtual void shutDown() = 0;
		virtual void init() = 0;
		virtual void simulate() = 0;
		virtual GLuint getRenderedTexture() = 0;
	};
}

#endif