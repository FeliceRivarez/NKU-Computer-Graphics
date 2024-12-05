#pragma once
#ifndef __INSPECTOR_VIEW_H__
#define __INSPECTOR_VIEW_H__

#include "glad/glad.h"
#include "glfw3.h"

#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Configure.h"
#include "Manager.h"
#include "Logger.h"
#include <iostream>
#include <string>

namespace FluidSimulation {

	class InspectorView {
	private:
		GLFWwindow* window;
		ImVec2 pos;

	public:

		int showID;

		InspectorView();
		InspectorView(GLFWwindow* window);
		void display();
	};
}

#endif