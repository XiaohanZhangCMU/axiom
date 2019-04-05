get_filename_component(target ${CMAKE_CURRENT_SOURCE_DIR} NAME)   # put current folder name as target name
include_directories(${include_directories} "include")
file (GLOB SRCS "src/*.cpp")
pybind11_add_module(${target} ${SRCS})
set_target_properties(${target} PROPERTIES PREFIX "" LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../lib/")
