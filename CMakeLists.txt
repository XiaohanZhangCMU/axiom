cmake_minimum_required(VERSION 3.4.0)
find_package(PythonLibs ${USE_PYTHON_VERSION} REQUIRED)
find_package(PythonInterp ${USE_PYTHON_VERSION} REQUIRED)

project (AXIOM)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR "\nIn-source builds are not allowed.  Clean up by:\n rm -rf CMakeCache.txt CMakeFiles")
endif()

site_name(HostName)
message("HostName = ${HostName}")

set(HOME $ENV{HOME})
message("HOME = ${HOME}")

set(PYTHON_LIBRARY "/Users/x/miniconda3/lib")
find_package(PythonLibs)

include_directories(${include_directories} "include" ${CPLUS_INCLUDE_DIR} ${PYTHON_INCLUDE_DIR})
add_subdirectory(pybind11)

# These lines are machine dependent.
if("${HostName}" MATCHES "rice*" )
    message("On rice: add self-defined blas and cuda library")
    include_directories(${include_directories} /usr/local/opt/openblas/include)
    link_directories   (${link_directories} /usr/local/opt/openblas/lib /usr/local/cuda/lib64)
endif()
if(APPLE)
    message("On OSX: add self-defined blas library")
    include_directories(${include_directories} /usr/local/opt/openblas/include)
    link_directories   (${link_directories} /usr/local/opt/openblas/lib)
endif()

# --------------------- CUDA  ---------------------- #

find_package(CUDA)
if (APPLE)
    add_definitions(-DCPU_ONLY)
    set(CPU_ONLY "On")
    message("Turn off CUDA on OSX.")
else()
    if(CUDA_FOUND)
        option (CUDA_ENFORCE_HOST_COMPILER "Force nvcc to use the same compiler used to compile .c(pp) files insted of gcc/g++" ON)
        if(${HostName} MATCHES "sh-*")
            set(CUDA_ARCH -gencode;arch=compute_61,code=sm_61)
            message("On Sherlock2.0. Use code=sm_61")
        else()
            set(CUDA_ARCH -gencode;arch=compute_35,code=sm_35)
            message("On Rice or Sherlock1.0. Use code=sm_35")
        endif()

        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_ARCH} -std=c++11")
        file (GLOB AXIOM_CUDA_SRCS "src/*.cu" "src/tensor.cpp")
        cuda_add_library(cuda_axiom SHARED ${AXIOM_CUDA_SRCS})

        target_link_libraries(cuda_axiom cublas)
        set_target_properties(cuda_axiom PROPERTIES LIBRARY_OUTPUT_DIRECTORY "../lib")
    else()
        add_definitions(-DCPU_ONLY)
        set(CPU_ONLY "On")
        message("Unable to find CUDA.")
    endif()
endif()

file (GLOB AXIOM_SRCS "src/*.cpp")
if(UNIX OR APPLE OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -std=c++11")
else()
    message(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.= ${CMAKE_COMPILER_IS_GNUCXX}")
endif()

#target_link_libraries(axiomlib openblas python3.6m)
#set_target_properties(axiomlib PROPERTIES PREFIX ""  LIBRARY_OUTPUT_DIRECTORY "../lib")

if (CUDA_FOUND)
    if(NOT CPU_ONLY)
        target_link_libraries(axiomlib cuda_axiom cublas )
    endif()
endif()

# compile shared library to interface with python
pybind11_add_module(axiomlib ${AXIOM_SRCS})


# -------------------- clean --------------------- #

add_custom_target(clean-all
    command ${CMAKE_MAKE_PROGRAM} clean
    command rm -rf CMakeFiles CMakeCache.txt cmake*.cmake
)

