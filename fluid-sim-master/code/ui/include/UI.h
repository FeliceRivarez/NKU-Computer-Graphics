#pragma once
#ifndef __UI_H__
#define __UI_H__

#include "glad/glad.h"
#include "glfw3.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Manager.h"
#include "Logger.h"

namespace FluidSimulation 
{
	using namespace std;

	class UI {
	private:
	public:
		UI();
		void run();
	};
}

#endif