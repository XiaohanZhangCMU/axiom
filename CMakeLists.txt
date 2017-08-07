cmake_minimum_required(VERSION 2.8)
project (axiom)
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR "\nIn-source builds are not allowed.  Clean up by:\n rm -rf CMakeCache.txt CMakeFiles")
endif()

# ---------------- DIRS & INCLUDE ------------------ #

set(BOOST_ROOT "/home/xzhang11/boost")
set(BOOST_INCLUDE_DIR "${BOOST_ROOT}/include")
set(BOOST_LIBRARIES "${BOOST_ROOT}/lib")
include_directories("include" "src" /home/xzhang11/bin/miniconda3/lib/python3.6/site-packages/numpy/core/include ${CPLUS_INCLUDE_DIR} ${BOOST_INCLUDE_DIR})

# This line is needed on rice.
# link_directories(/usr/local/cuda/lib64)

# --------------------- CUDA  ---------------------- #

find_package(CUDA)
if(CUDA_FOUND)
    option (CUDA_ENFORCE_HOST_COMPILER "Force nvcc to use the same compiler used to compile .c(pp) files insted of gcc/g++" ON)
    set(CUDA_ARCH -gencode;arch=compute_35,code=sm_35)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_ARCH} -std=c++11")
    file (GLOB AXIOM_CUDA_SRCS "src/*.cu")
    message("axiom_cuda_src = ${AXIOM_CUDA_SRCS}, cuda_libr = ${CUDA_LIBRARIES}")
    cuda_add_library(cuda_axiom SHARED ${AXIOM_CUDA_SRCS})
    target_link_libraries(cuda_axiom ${Boost_LIBRARIES})
    set_target_properties(cuda_axiom PROPERTIES LIBRARY_OUTPUT_DIRECTORY "../lib")
elseif(not CUDA_FOUND)
    message(FATAL_ERROR "Unable to find CUDA. Did you set CUDA_ROOT?")
endif()

# ----------------- BOOST+axiom  ------------------ #

find_package(Boost 1.64.0 REQUIRED)
if(Boost_FOUND)
    set(Boost_USE_STATIC_LIBS OFF) 
    set(Boost_USE_MULTITHREADED ON)  
    set(Boost_USE_STATIC_RUNTIME OFF) 
    set(Boost_USE_STATIC_LIBS OFF)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_USE_STATIC_RUNTIME OFF)
    find_package(Boost 1.64.0 COMPONENTS python3 thread)
    file (GLOB AXIOM_SRCS "src/*.cpp")
    if(UNIX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
    else()
        message(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.")
    endif()
    add_library(axiom SHARED ${AXIOM_SRCS})
    message("Boost libr = ${Boost_LIBRARIES}.")
  
    target_link_libraries(axiom cuda_axiom ${Boost_LIBRARIES} cublas blas)
    set_target_properties(axiom PROPERTIES PREFIX ""  LIBRARY_OUTPUT_DIRECTORY "../lib")
    if(UNIX AND (CUDA_FOUND OR CMAKE_COMPILER_IS_GNUCXX))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -std=c++11")
    else()
        message(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.= ${CMAKE_COMPILER_IS_GNUCXX}")
    endif()
elseif(not Boost_FOUND)
    message(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
endif()

# -------------------- clean --------------------- #

add_custom_target(clean-all
    command ${CMAKE_MAKE_PROGRAM} clean
    command rm -rf CMakeFiles CMakeCache.txt cmake*.cmake
)

