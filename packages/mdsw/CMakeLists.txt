get_filename_component(target ${CMAKE_CURRENT_SOURCE_DIR} NAME)   # put current folder name as target name

if(APPLE)
    if(${HostName} MATCHES "Mohans-iMac.local") #imac
        set(X11_DIR "/usr/X11/")
    elseif(${HostName} MATCHES "mmc") #macbook
        set(X11_DIR "/usr/X11R6")
    else()
        set(X11_DIR "/opt/X11/")
    endif()
endif()
include_directories(${include_directories} "include" "../../ext_pybind11/include/" ${X11_include_dir} "${X11_DIR}/include")
link_directories(${link_directories} "${X11_DIR}/lib")

#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()

add_definitions(-DNO_XPM -D_SW_Si)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix -Wno-unused-but-set-variable -Wno-misleading-indentation -Wno-format-overflow -Wno-unused-variable -Wno-parentheses -Wno-restrict -Wno-class-memaccess")
file (GLOB SRCS "src/*.cpp")

set(EXECUTABLE OFF)
if(EXECUTABLE)
    add_executable(${target} ${SRCS})
    target_compile_options(${target} PUBLIC ${PYBIND11_CPP_STANDARD})
    set_target_properties(${target} PROPERTIES PREFIX "" RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/")
else()
    pybind11_add_module(${target} ${SRCS})
    set(OUTPUTDIR "${CMAKE_CURRENT_SOURCE_DIR}/../../lib")
    set_target_properties(${target} PROPERTIES PREFIX "" LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../lib/")
endif()
target_link_libraries(${target} PRIVATE X11)
