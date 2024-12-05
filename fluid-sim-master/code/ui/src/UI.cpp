#include "UI.h"

namespace FluidSimulation {

    UI::UI() {

    }

	void UI::run() {

        // create and initialize imgui
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGui::StyleColorsDark();
        ImGuiIO& io = ImGui::GetIO();
        io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
        (void)io;
        io.Fonts->AddFontFromFileTTF("C:\\Windows\\Fonts\\consola.ttf", fontSize);
        

        // initialize glfw
        if (!glfwInit())
        {
            Glb::Logger::getInstance().addLog("GLFW initialization failed.");
            return;
        }
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        // create window
        GLFWwindow* window = glfwCreateWindow(windowWidth, windowHeight, "Fluid Simulation System", NULL, NULL);
        if (!window) {
            Glb::Logger::getInstance().addLog("Fail to create window.");
            glfwTerminate();
            return;
        }
        // set it as current context
        glfwMakeContextCurrent(window);

        // Maximize
        // glfwMaximizeWindow(window);


        // initialize glad
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
            Glb::Logger::getInstance().addLog("GLAD initialization failed.");
            return;
        }
        Glb::Logger::getInstance().addLog("GLFW and GLAD initialization succeeded.");

        // make imgui ready for glfw and OpenGL
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 330");
        ImGui_ImplOpenGL3_CreateDeviceObjects();

        // init manager
        Manager::getInstance().init(window);
        Glb::Logger::getInstance().addLog("Start Main Render Loop...");

        // main render loop
        while (!glfwWindowShouldClose(window)) {
            
            // deal with events
            glfwPollEvents();

            // 
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            
            // render
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            Manager::getInstance().displayViews();

            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

            // refresh
            glfwSwapBuffers(window);
        }

        // release
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
	}
}