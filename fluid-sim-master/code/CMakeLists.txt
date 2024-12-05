
cmake_minimum_required (VERSION 3.8)

project(FluidSimulationSystem VERSION 1.0
		DESCRIPTION "A simple system for fluid simulation"
		LANGUAGES C CXX)

# where to find the .h
include_directories(
	"./third_party/imgui/include"
	"./third_party/glfw/include"
	"./third_party/glad/include"
	"./third_party/boost_1_84_0"
	"./third_party/glm"
	"./third_party/stb"
	"./common/include"
	"./fluid2d/"
	"./fluid2d/Lagrangian/include"
	"./fluid2d/Eulerian/include"
	"./fluid3d/"
	"./fluid3d/Lagrangian/include"
	"./fluid3d/Eulerian/include"
	"./ui/include"
	"."
)

# third party
add_subdirectory("./third_party/imgui")
add_subdirectory("./third_party/glad")
add_subdirectory("./third_party/glm")

# common
add_subdirectory("./common")

# 2d simulation scheme
add_subdirectory("./fluid2d/Lagrangian")
add_subdirectory("./fluid2d/Eulerian")

# 3d simulation scheme
add_subdirectory("./fluid3d/Lagrangian")
add_subdirectory("./fluid3d/Eulerian")

# ui
add_subdirectory("./ui")

# exe
add_executable (FluidSimulationSystem "code.cpp" "code.h")

# just need ui's lib to run
target_link_libraries(FluidSimulationSystem PRIVATE ui)
# target_include_directories(${PROJECT_NAME} PRIVATE "./ui/include")


