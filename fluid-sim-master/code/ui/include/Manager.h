#pragma once
#ifndef __MANAGER_H__
#define __MANAGER_H__

#include "glad/glad.h"
#include "glfw3.h"

#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "SceneView.h"
#include "InspectorView.h"
#include "ProjectView.h"

#include "Lagrangian2dComponent.h"
#include "Eulerian2dComponent.h"
#include "Lagrangian3dComponent.h"
#include "Eulerian3dComponent.h"

#include <vector>

namespace FluidSimulation
{
	class SceneView;
	class InspectorView;
	class ProjectView;

	class Manager {
	public:
		static Manager& getInstance() {
			static Manager instance;
			return instance;
		}

		void init(GLFWwindow* window);
		void displayViews();
		void displayToolBar();

		SceneView* getSceneView() const { return sceneView; };
		InspectorView* getInspectorView() const { return inspectorView; };
		ProjectView* getProjectView() const { return projectView; };
		GLFWwindow* getWindow() const { return window; };
		Glb::Component* getMethod() const { return currentMethod; };
		void setMethod(Glb::Component* method) { currentMethod = method; };

	private:

		Manager() {
			window = NULL;
			sceneView = NULL;
			inspectorView = NULL;
			projectView = NULL;
			currentMethod = NULL;
		};

		Manager(const Manager&) = delete;
		Manager& operator=(const Manager&) = delete;
		
		GLFWwindow* window;
		SceneView* sceneView;
		InspectorView* inspectorView;
		ProjectView* projectView;

		Glb::Component* currentMethod;
		
	};
}

#endif