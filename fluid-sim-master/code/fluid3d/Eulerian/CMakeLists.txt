
enable_language(C CXX)

file(GLOB_RECURSE Eulerian3D_SOURCE_FILES "./src/*.cpp")
file(GLOB_RECURSE Eulerian3D_HEADER_FILES "./include/*.h ./include/*.hpp")

source_group("Header Files" FILES ${Eulerian3D_HEADER_FILES})

add_library(eulerian3d STATIC "${Eulerian3D_SOURCE_FILES}" "${Eulerian3D_HEADER_FILES}")
target_include_directories(eulerian3d PRIVATE "./include")

include_directories("./include")

# common
target_link_libraries(eulerian3d common)

# glfw
target_link_libraries(eulerian3d "${PROJECT_SOURCE_DIR}/third_party/glfw/lib/glfw3.lib")