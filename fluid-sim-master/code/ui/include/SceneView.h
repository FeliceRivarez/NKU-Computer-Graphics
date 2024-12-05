#pragma once
#ifndef __SCENE_VIEW_H__
#define __SCENE_VIEW_H__

#include <iostream>

#include "glad/glad.h"
#include "glfw3.h"

#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Configure.h"
#include "Component.h"

#include "Global.h"

#include "Manager.h"
#include <thread>


namespace FluidSimulation {
	class SceneView {

	public:
		GLFWwindow* window;
		ImVec2 pos;

		// std::thread renderThread;

		GLuint texture;

		// 检测鼠标是否停留在渲染界面
		bool isMouseHovering = false;
		// 鼠标拖动状态
		bool isLeftMouseDragging = false;
		bool isRightMouseDragging = false;
		double lastMouseX = 0.0;
		double lastMouseY = 0.0;
		double mouseX, mouseY;

	public:
		SceneView();
		SceneView(GLFWwindow* window);
		void display();
	};
}

#endif