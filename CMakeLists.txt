cmake_minimum_required(VERSION 3.4.0)
site_name(HostName)
message("HostName = ${HostName}")
set(HOME $ENV{HOME})
message("HOME = ${HOME}")

set(CMAKE_C_COMPILER /usr/local/bin/gcc)
set(CMAKE_CXX_COMPILER /usr/local/bin/g++)
SET(CMAKE_C_LINK_EXECUTABLE /usr/bin/ld)


project (AXIOM)


if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR "\nIn-source builds are not allowed.  Clean up by:\n rm -rf CMakeCache.txt CMakeFiles")
endif()

set(PYTHON_LIBRARY "/Users/x/miniconda3/lib")
find_package(PythonLibs 3.6 REQUIRED)
find_package(PythonInterp 3.6 REQUIRED)

set(BOOST_ROOT "${HOME}/boost")
set(BOOST_INCLUDE_DIR "${BOOST_ROOT}/include")
set(BOOST_LIBRARIES_DIR "${BOOST_ROOT}/lib")
find_package(Boost 1.68.0 COMPONENTS python36 thread system REQUIRED)

# --------------- include and link directories  ---------------- #

include_directories(${include_directories} ${CPLUS_INCLUDE_DIR} ${PYTHON_INCLUDE_DIR})
add_subdirectory(pybind11)

# These lines are machine dependent.
if("${HostName}" MATCHES "rice*" )
    message("On rice: add self-defined blas and cuda library")
    include_directories(${include_directories} /usr/local/opt/openblas/include)
    link_directories   (${link_directories} ${PYTHON_LIBRARY} /usr/local/opt/openblas/lib /usr/local/cuda/lib64 ${BOOST_LIBRARIES_DIR})
endif()
if(APPLE)
    # message("On OSX: add self-defined blas library")
    include_directories(${include_directories} /usr/local/opt/openblas/include)
    link_directories   (${link_directories} ${PYTHON_LIBRARY} /usr/local/opt/openblas/lib ${BOOST_LIBRARIES_DIR} )
endif()



# This part also needs to be moved to each module folder, because some module may not use CUDA
# --------------------------  CUDA ---------------------------- #

find_package(CUDA)
if (APPLE)
    # message("Turn off CUDA on OSX.")
    add_definitions(-DCPU_ONLY)
    set(CPU_ONLY "On")
else()
    if(CUDA_FOUND)
        option (CUDA_ENFORCE_HOST_COMPILER "Force nvcc to use the same compiler used to compile .c(pp) files insted of gcc/g++" ON)
        if(${HostName} MATCHES "sh-*")
            # message("On Sherlock2.0. Use code=sm_61")
            set(CUDA_ARCH -gencode;arch=compute_61,code=sm_61)
        else()
            # message("On Rice or Sherlock1.0. Use code=sm_35")
            set(CUDA_ARCH -gencode;arch=compute_35,code=sm_35)
        endif()

        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_ARCH} -std=c++11")
        file (GLOB AXIOM_CUDA_SRCS "src/*.cu" "src/tensor.cpp")
        cuda_add_library(cuda_axiom SHARED ${AXIOM_CUDA_SRCS})

        target_link_libraries(cuda_axiom cublas)
        set_target_properties(cuda_axiom PROPERTIES LIBRARY_OUTPUT_DIRECTORY "../lib")
    else()
        # message("Unable to find CUDA.")
        add_definitions(-DCPU_ONLY)
        set(CPU_ONLY "On")
    endif()
endif()

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
ADD_SUBDIRECTORY(zoo)
ADD_SUBDIRECTORY(tensor)


# -------------------- clean --------------------- #

add_custom_target(clean-all
    command ${CMAKE_MAKE_PROGRAM} clean
    command rm -rf CMakeFiles CMakeCache.txt cmake*.cmake
)

