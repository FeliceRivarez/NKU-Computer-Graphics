#include "SceneView.h"

namespace FluidSimulation {

	SceneView::SceneView() {

	}

	SceneView::SceneView(GLFWwindow* window) {
		this->window = window;
        glfwMakeContextCurrent(window);
        texture = -1;   // no texture;
	}

	void SceneView::display() {

		ImGui::Begin("Scene", NULL, ImGuiWindowFlags_NoCollapse);

        ImVec2 windowSize = ImGui::GetWindowSize();

        if (Manager::getInstance().getMethod() != NULL) {
            if (simulating || texture == -1) {
                Manager::getInstance().getMethod()->simulate();
                texture = Manager::getInstance().getMethod()->getRenderedTexture();
                Glb::Timer::getInstance().timeFPS();
            }
            else if (!simulating) {
                texture = Manager::getInstance().getMethod()->getRenderedTexture();
                Glb::Timer::getInstance().timeFPS();
            }
            ImGui::Text(("FPS: " + Glb::Timer::getInstance().getFPS()).c_str());
        }


        // show image
        ImGui::SetCursorPosX((windowSize.x - imageWidth) * 0.5f);
        ImGui::SetCursorPosY(55);
        ImVec2 imageSize(imageWidth, imageHeight);
        if (Manager::getInstance().getMethod() == NULL) {
            ImGui::Image(NULL, imageSize);
        }
        else { 
            ImGui::Image((ImTextureID)texture, imageSize); 
        }
        ImGui::GetWindowDrawList()->AddImage((ImTextureID)texture,
            ImVec2(ImGui::GetItemRectMin().x, ImGui::GetItemRectMin().y), ImVec2(ImGui::GetItemRectMax().x, ImGui::GetItemRectMax().y),
            ImVec2(0, 1), ImVec2(1, 0));

        if (ImGui::IsItemHovered() && Manager::getInstance().getMethod() != NULL) {

            glfwGetCursorPos(window, &mouseX, &mouseY);

            Glb::Camera::getInstance().ProcessScale(static_cast<float>(ImGui::GetIO().MouseWheel));

            if (ImGui::IsMouseDragging(0, 0.0f)) {
                if (!isLeftMouseDragging) {
                    isLeftMouseDragging = true;
                    lastMouseX = mouseX;
                    lastMouseY = mouseY;
                }
                float deltaX = (float)(mouseX - lastMouseX);
                float deltaY = (float)(mouseY - lastMouseY);

                Glb::Camera::getInstance().ProcessRotate(glm::vec2(deltaX, deltaY));
                
                lastMouseX = mouseX;
                lastMouseY = mouseY;
            }
            else {
                isLeftMouseDragging = false;
            }

            if (ImGui::IsMouseDragging(1, 0.0f)) {
                if (!isRightMouseDragging) {
                    isRightMouseDragging = true;
                    lastMouseX = mouseX;
                    lastMouseY = mouseY;
                }
                float deltaX = (float)(mouseX - lastMouseX);
                float deltaY = (float)(mouseY - lastMouseY);

                Glb::Camera::getInstance().ProcessMove(glm::vec2(deltaX, deltaY));

                lastMouseX = mouseX;
                lastMouseY = mouseY;
            }
            else {
                isRightMouseDragging = false;
            }
        }

		ImGui::End();
	}
}