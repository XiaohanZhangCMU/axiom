# put current folder name as target name
get_filename_component(target ${CMAKE_CURRENT_SOURCE_DIR} NAME) 
set(cuda_target "cuda_${target}")

set(CPU_ONLY "Off")

if (APPLE) 
    # message("Always turn off CUDA on OSX.")
    set(CPU_ONLY "On")
endif()
if (CPU_ONLY) 
    add_definitions(-DCPU_ONLY)
    message("Add defintion cpu_only")
endif()

set(BOOST_ROOT "${HOME}/boost")
set(BOOST_INCLUDE_DIR "${BOOST_ROOT}/include")
set(BOOST_LIBRARIES_DIR "${BOOST_ROOT}/lib")
message("BOOST_ROOT = ${BOOST_ROOT}")

# These lines are machine dependent.
set(BLAS_INCLUDE_DIR "/usr/local/opt/openblas/include")
set(BLAS_LIB_DIR "/usr/local/opt/openblas/lib")
set(CUBLAS_LIB_DIR "/usr/local/cuda/lib64")

if(APPLE)
    find_package(Boost 1.66.0 COMPONENTS python36 thread system REQUIRED)
else()
    find_package(Boost 1.66.0 COMPONENTS python3 thread system REQUIRED)
endif()

include_directories(${include_directories} ${BLAS_INCLUDE_DIR} ${BOOST_INCLUDE_DIR} "include" "../pybind11/include/")
link_directories   (${link_directories} ${PYTHON_LIBRARY} ${CUBLAS_LIB_DIR} ${BLAS_LIB_DIR} ${BOOST_LIBRARIES_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/../lib")

#set(CUDA_ROOT "/share/software/user/open/cuda/9.0.176")
if(NOT CPU_ONLY)
    find_package(CUDA REQUIRED)
    message("I am here 1")
    enable_language("CUDA")
    option (CUDA_ENFORCE_HOST_COMPILER "Force nvcc to use the same compiler used to compile .c(pp) files insted of gcc/g++" ON)
    if(${HostName} MATCHES "sh-*")
        # message("On Sherlock2.0. Use code=sm_61")
        set(CUDA_ARCH -gencode;arch=compute_61,code=sm_61)
    else()
        # message("On Rice or Sherlock1.0. Use code=sm_35")
        set(CUDA_ARCH -gencode;arch=compute_35,code=sm_35)
    endif()

    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_ARCH} -std=c++11")
     
    file (GLOB AXIOM_CUDA_SRCS "src/*.cu")
    #cuda_add_library(cuda_axiom SHARED ${AXIOM_CUDA_SRCS})
    add_library(${cuda_target} SHARED ${AXIOM_CUDA_SRCS})
    #add_executable(${cuda_target} ${AXIOM_CUDA_SRCS})
    set_target_properties(${cuda_target} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(${cuda_target} cublas ${CUDA_LIBRARIES})
    set_target_properties(${cuda_target} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../lib")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix -Wno-unused-but-set-variable -Wno-misleading-indentation -Wno-format-overflow -Wno-unused-variable -Wno-parentheses -Wno-restrict -Wno-class-memaccess")
file (GLOB SRCS "src/*.cpp")
message("SRCS = ${SRCS}")
pybind11_add_module(${target} ${SRCS})
#pybind11_add_exe(${target} ${SRCS})
message("cuda_libraries = ${CUDA_LIBRARIES}")
if(NOT CPU_ONLY)
    target_link_libraries(${target} PRIVATE openblas ${cuda_target} ${CUDA_LIBRARIES})
else()
    target_link_libraries(${target} PRIVATE openblas)
endif()

set_target_properties(${target} PROPERTIES PREFIX "" LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../lib/")


