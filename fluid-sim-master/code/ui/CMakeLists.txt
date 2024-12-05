# ui/CMakeLists.txt
enable_language(C CXX)

file(GLOB_RECURSE UI_SOURCE_FILES "./src/*.cpp")
file(GLOB_RECURSE UI_HEADER_FILES "./include/*.h ./include/*.hpp")

source_group("Header Files" FILES ${UI_HEADER_FILES})

add_library(ui STATIC "${UI_SOURCE_FILES}" "${UI_HEADER_FILES}")
target_include_directories(ui PRIVATE "./include")

# glfw
target_link_libraries(ui glad)
target_link_libraries(ui imgui)
target_link_libraries(ui lagrangian2d)
target_link_libraries(ui eulerian2d)
target_link_libraries(ui lagrangian3d)
target_link_libraries(ui eulerian3d)

target_link_libraries(ui "${PROJECT_SOURCE_DIR}/third_party/glfw/lib/glfw3.lib")