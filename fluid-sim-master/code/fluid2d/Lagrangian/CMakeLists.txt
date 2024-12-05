
enable_language(C CXX)

file(GLOB_RECURSE Lagrangian2D_SOURCE_FILES "./src/*.cpp")
file(GLOB_RECURSE Lagrangian2D_HEADER_FILES "./include/*.h ./include/*.hpp")

source_group("Header Files" FILES ${Lagrangian2D_HEADER_FILES})

add_library(lagrangian2d STATIC "${Lagrangian2D_SOURCE_FILES}" "${Lagrangian2D_HEADER_FILES}")
target_include_directories(lagrangian2d PRIVATE "./include")

# common
target_link_libraries(lagrangian2d common)

# glfw
target_link_libraries(lagrangian2d "${PROJECT_SOURCE_DIR}/third_party/glfw/lib/glfw3.lib")