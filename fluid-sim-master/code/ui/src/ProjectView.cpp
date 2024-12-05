#include "ProjectView.h"

namespace FluidSimulation {
	ProjectView::ProjectView() {

	}

	ProjectView::ProjectView(GLFWwindow* window) {
		this->window = window;

	}

	void ProjectView::display() {
		ImGui::Begin("Project", NULL, ImGuiWindowFlags_NoCollapse);

		ImGui::BeginChild("ScrollingRegion", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar);

		for (const auto& message : Glb::Logger::getInstance().getLog()) {
			ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));
			ImGui::TextWrapped("%s", message.c_str());
			ImGui::PopStyleColor();
		}

		if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
			ImGui::SetScrollHereY(1.0f);

		ImGui::EndChild();
		ImGui::End();
	}
}