#pragma once
#ifndef __PROJECT_VIEW_H__
#define __PROJECT_VIEW_H__

#include "glad/glad.h"
#include "glfw3.h"

#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Manager.h"

namespace FluidSimulation {
	class ProjectView {
	private:
		GLFWwindow* window;
		ImVec2 pos;
	public:
		ProjectView();
		ProjectView(GLFWwindow* window);
		void display();
	};
}

#endif