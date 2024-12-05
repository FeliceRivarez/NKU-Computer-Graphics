enable_language(C CXX)

file(GLOB_RECURSE COMMON_SOURCE_FILES "./src/*.cpp")
file(GLOB_RECURSE COMMON_HEADER_FILES "./include/*.h ./include/*.hpp")

source_group("Header Files" FILES ${COMMON_HEADER_FILES})

add_library(common STATIC "${COMMON_SOURCE_FILES}" "${COMMON_HEADER_FILES}")
target_include_directories(common PRIVATE "./include")