get_filename_component(target ${CMAKE_CURRENT_SOURCE_DIR} NAME)   # put current folder name as target name
include_directories(${include_directories} "include" "../pybind11/include/")
message("${include_directories}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix -Wno-unused-but-set-variable -Wno-misleading-indentation -Wno-format-overflow -Wno-unused-variable -Wno-parentheses -Wno-restrict -Wno-class-memaccess")
file (GLOB SRCS "src/*.cpp")
pybind11_add_module(${target} ${SRCS})
target_link_libraries(${target} PRIVATE openblas)
set_target_properties(${target} PROPERTIES PREFIX "" LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../lib/")





# set(target "tmpexe")
# #add_executable(${target} ${SRCS})
# add_executable(${target} src/main.cpp ${SRCS})
# target_compile_options(${target} PUBLIC ${PYBIND11_CPP_STANDARD})
# #pybind11_add_exe(${target} ${SRCS})
# # add_definitions(-c)
#
# #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -c")
#
# #set_target_properties(${target} PROPERTIES PREFIX "" RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin/")
# #set_target_properties(${target} PROPERTIES PREFIX "" RUNTIME_OUTPUT_DIRECTORY "./")
# set(EXECUTABLE_OUTPUT_PATH /Users/x/)
#
#
# message("runtime=${RUNTIME_OUTPUT_DIRECTORY}")
# message("boost = ${Boost_LIBRARIES}")
# message("curr = ${CMAKE_CURRENT_SOURCE_DIR}")
# target_link_libraries(${target} PRIVATE openblas python3.6m ${Boost_LIBRARIES})
#
# # # compile executable
# # # pybind11_add_exe(${target}.exe ${SRCS})
# # add_executable(${target}.exe ${SRC})
# # target_include_directories(${target_name}.exe "${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include")
# # target_link_libraries(${target}.exe PRIVATE openblas python3.6m ${Boost_LIBRARIES})
# # set_target_properties(${target}.exe PROPERTIES PREFIX ""  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin/")
# #
# # # libc++ or libstdc++? That is a question to resolve
# # SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++")

