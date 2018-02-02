cmake_minimum_required(VERSION 3.4.0)
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
set(BOOST_LIBRARIES_DIR "${BOOST_ROOT}/lib")

set(CONDA_ROOT ${HOME}/miniconda3)
set(CONDA_INCLUDE_DIR "${CONDA_ROOT}/include")
set(CONDA_LIBRARIES_DIR "${CONDA_ROOT}/lib")

include_directories(${include_directories} "include" ${CPLUS_INCLUDE_DIR} ${BOOST_INCLUDE_DIR})
link_directories   (${link_directories}    "lib" ${CONDA_LIBRARIES_DIR} ${BOOST_LIBRARIES_DIR})

# These lines are machine dependent.
if("${HostName}" MATCHES "rice*")
    message("On rice: add self-defined blas and cuda library")
    include_directories(${include_directories} /usr/local/opt/openblas/include)
    link_directories   (${link_directories} /usr/local/opt/openblas/lib /usr/local/cuda/lib64)
endif()

# --------------------- CUDA  ---------------------- #

if (APPLE)
    add_definitions(-DCPU_ONLY)
    set(CPU_ONLY "On")
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
            message("On Rice or Sherlock1.0. Use code=sm_35")
        endif()
    
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_ARCH} -std=c++11")
        file (GLOB AXIOM_CUDA_SRCS "src/*.cu" "src/tensor.cpp")
        cuda_add_library(cuda_axiom SHARED ${AXIOM_CUDA_SRCS})

        target_link_libraries(cuda_axiom cublas)
        set_target_properties(cuda_axiom PROPERTIES LIBRARY_OUTPUT_DIRECTORY "../lib")
    else() # if(not CUDA_FOUND)
        add_definitions(-DCPU_ONLY)
        set(CPU_ONLY "On")
        message("Unable to find CUDA. Set CPU_ONLY.")
    endif()
endif()

# ----------------- BOOST+axiom  ------------------ #

find_package(Boost 1.66.0 REQUIRED)
if(Boost_FOUND)
    set(Boost_USE_STATIC_LIBS OFF) 
    set(Boost_USE_MULTITHREADED ON)  
    set(Boost_USE_STATIC_RUNTIME OFF)
    #find_package(Boost 1.66.0 COMPONENTS python3 thread system numpy3)
    find_package(Boost 1.66.0 COMPONENTS python3 thread system)

    file (GLOB AXIOM_SRCS "src/*.cpp")
    if(UNIX OR APPLE OR CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -std=c++11")
    else()
        message(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.= ${CMAKE_COMPILER_IS_GNUCXX}")
    endif()

    # compile shared library to interface with python
    add_library(axiomlib SHARED ${AXIOM_SRCS})
    if (NOT CPU_ONLY)
        target_link_libraries(axiomlib cuda_axiom cublas )
    endif()
    target_link_libraries(axiomlib openblas python3.6m ${Boost_LIBRARIES})
    set_target_properties(axiomlib PROPERTIES PREFIX ""  LIBRARY_OUTPUT_DIRECTORY "../lib")

    # compile executable
    add_executable(axiomexe ${AXIOM_SRCS})
    if (NOT CPU_ONLY)
        target_link_libraries(axiomexe cuda_axiom cublas )
    endif()
    target_link_libraries(axiomexe openblas python3.6m ${Boost_LIBRARIES})
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

