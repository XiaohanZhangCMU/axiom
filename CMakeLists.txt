# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# MAIN CMAKElISTS.TXT TO COMPILE SUBDIRECTORIES IN ROOT FOLDER OF AXIOM
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.1.0)
site_name(HostName)
message("HostName = ${HostName}")
set(HOME $ENV{HOME})
message("HOME = ${HOME}")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# CHANGE THE FOLLOWING SETTINGS OR COMMENT OUT TO USE DEFAULTS.
#           BETTER DEFINE DIFFERENT MACHINES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if(APPLE)
    set(CMAKE_C_COMPILER /usr/local/bin/gcc)
    set(CMAKE_CXX_COMPILER /usr/local/bin/g++)
    SET(CMAKE_C_LINK_EXECUTABLE /usr/bin/ld)
    set(PYTHON_LIBRARY "/Users/x/miniconda3/lib")
else()
    #set(CMAKE_C_COMPILER /share/software/user/open/gcc/8.1.0/bin/gcc)
    #set(CMAKE_CXX_COMPILER /share/software/user/open/gcc/8.1.0/bin/g++)
    #SET(CMAKE_C_LINK_EXECUTABLE /usr/bin/ld)
    set(PYTHON_LIBRARY "/scratch/users/xzhang11/miniconda3/lib")
    set(PYTHON_INCLUDE_DIR "/scratch/users/xzhang11/miniconda3/include/python3.6m")
endif()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# FIND PYTHON LIBRARIES AND COMPILE PYBIND11
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

project (AXIOM)
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR "\nIn-source builds are not allowed.  Clean up by:\n rm -rf CMakeCache.txt CMakeFiles")
endif()

find_package(PythonLibs 3.6 REQUIRED)
find_package(PythonInterp 3.6 REQUIRED)
include_directories(${include_directories} ${CPLUS_INCLUDE_DIR} ${PYTHON_INCLUDE_DIR})
add_subdirectory(pybind11)
if(UNIX OR APPLE OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -std=c++11")
else()
    message(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.= ${CMAKE_COMPILER_IS_GNUCXX}")
endif()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# EXAMPLES (ADD YOUR OWN DIRECTORIES BELOW)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

ADD_SUBDIRECTORY(zoo)
ADD_SUBDIRECTORY(tensor)
ADD_SUBDIRECTORY(mdfem)
ADD_SUBDIRECTORY(mdsw)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# CLEAN
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

add_custom_target(clean-all
    command ${CMAKE_MAKE_PROGRAM} clean
    command rm -rf CMakeFiles CMakeCache.txt cmake*.cmake
)

