cmake_minimum_required(VERSION 2.8)
project (AXIOM)
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR "\nIn-source builds are not allowed.  Clean up by:\n rm -rf CMakeCache.txt CMakeFiles")
endif()

# ---------------- DIRS & INCLUDE ------------------ #
site_name(HostName)
message("HostName = ${HostName}")
set(HOME $ENV{HOME})
message("HOME = ${HOME}")
set(BOOST_ROOT "${HOME}/boost")
set(BOOST_INCLUDE_DIR "${BOOST_ROOT}/include")
set(BOOST_LIBRARIES "${BOOST_ROOT}/lib")
include_directories("include" "src" ${HOME}/miniconda3/lib/python3.6/site-packages/numpy/core/include /usr/local/opt/openblas/include ${CPLUS_INCLUDE_DIR} ${BOOST_INCLUDE_DIR})

# This line is needed on rice to find cublas.
link_directories(${link_directories} /usr/local/opt/openblas/lib ${HOME}/miniconda3/lib)

# --------------------- CUDA  ---------------------- #

if (APPLE)
    add_definitions(-DCPU_ONLY)
    message("Always set CPU_ONLY on OSX. Turn this off if your computer has GPU cards.")
else()
    find_package(CUDA)
    if(CUDA_FOUND)
        option (CUDA_ENFORCE_HOST_COMPILER "Force nvcc to use the same compiler used to compile .c(pp) files insted of gcc/g++" ON)
    
        if(${HostName} MATCHES "sh-*")
            set(CUDA_ARCH -gencode;arch=compute_61,code=sm_61)
            message("On Sherlock2.0. Use code=sm_61")
        else()
            set(CUDA_ARCH -gencode;arch=compute_35,code=sm_35)
            message("On Rice/Sherlock1.0. Use code=sm_35")
        endif()
    
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_ARCH} -std=c++11")
        file (GLOB AXIOM_CUDA_SRCS "src/*.cu" "src/tensor.cpp")
        cuda_add_library(cuda_axiom SHARED ${AXIOM_CUDA_SRCS})
        link_directories( ${link_directories} /usr/local/cuda/lib64  )
        target_link_libraries(cuda_axiom ${Boost_LIBRARIES} cublas)
        set_target_properties(cuda_axiom PROPERTIES LIBRARY_OUTPUT_DIRECTORY "../lib")
    else() # if(not CUDA_FOUND)
        add_definitions(-DCPU_ONLY)
        message("Unable to find CUDA. Set CPU_ONLY.")
    endif()
endif()

# ----------------- BOOST+axiom  ------------------ #

find_package(Boost 1.66.0 REQUIRED)
if(Boost_FOUND)
    set(Boost_USE_STATIC_LIBS OFF) 
    set(Boost_USE_MULTITHREADED ON)  
    set(Boost_USE_STATIC_RUNTIME OFF)
    find_package(Boost 1.66.0 COMPONENTS python3 thread system numpy3)

    file (GLOB AXIOM_SRCS "src/*.cpp")
    if(UNIX OR APPLE OR CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -std=c++11")
    else()
        message(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.= ${CMAKE_COMPILER_IS_GNUCXX}")
    endif()

    add_library(axiomlib SHARED ${AXIOM_SRCS})

    if (CUDA_FOUND)
        target_link_libraries(axiomlib cuda_axiom ${Boost_LIBRARIES} cublas )
    endif()
    target_link_libraries(axiomlib ${Boost_LIBRARIES} openblas python3.6m)
    set_target_properties(axiomlib PROPERTIES PREFIX ""  LIBRARY_OUTPUT_DIRECTORY "../lib")

    add_executable(axiomexe ${AXIOM_SRCS})
    target_link_libraries(axiomexe ${Boost_LIBRARIES} openblas python3.6m)
    set_target_properties(axiomexe PROPERTIES PREFIX ""  RUNTIME_OUTPUT_DIRECTORY "../bin")

    # libc++ or libstdc++? That is a question to resolve
    # SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++")
elseif(not Boost_FOUND)
    message(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
endif()

# -------------------- clean --------------------- #

add_custom_target(clean-all
    command ${CMAKE_MAKE_PROGRAM} clean
    command rm -rf CMakeFiles CMakeCache.txt cmake*.cmake
)

