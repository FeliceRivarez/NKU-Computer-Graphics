
enable_language(C CXX)

file(GLOB_RECURSE Eulerian2D_SOURCE_FILES "./src/*.cpp")
file(GLOB_RECURSE Eulerian2D_HEADER_FILES "./include/*.h ./include/*.hpp")

source_group("Header Files" FILES ${Eulerian2D_HEADER_FILES})

add_library(eulerian2d STATIC "${Eulerian2D_SOURCE_FILES}" "${Eulerian2D_HEADER_FILES}")
target_include_directories(eulerian2d PRIVATE "./include")

include_directories("./include")

# common
target_link_libraries(eulerian2d common)

# glfw
target_link_libraries(eulerian2d "${PROJECT_SOURCE_DIR}/third_party/glfw/lib/glfw3.lib")