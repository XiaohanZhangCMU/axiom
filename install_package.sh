#!/bin/bash

RED='\033[0;31m'
NOCOLOR='\033[0m'

unameOut="$(uname -s)"

case ${unameOut} in
    Linux*) machine=Linux;;
    Darwin*) machine=Mac;;
    CYGWIN*) machine=Cygwin;;
    MINGW*) machine="UNKNOWN::${unameOut}"
esac
echo "Install package to using ${machine} configurations"

if [ ! -d build ];
then
    mkdir build
fi
cd build

cp ../CMakeLists.template ./CMakeLists.txt

valid_packages=$(echo $(ls ../packages) | tr " " "\n")

if [ -z $1 ]; then
    echo "${RED}Add all packages?${NOCOLOR}"
    for pkg in $valid_packages
    do
        echo "ADD_SUBDIRECTORY(packages/${pkg})" >> CMakeLists.txt
    done
    mv CMakeLists.txt ..
else
    Any="False"
    for request in "$@"
    do
        One="False"
        for pkg in $valid_packages
        do
            if [ $request == $pkg ] ; then
                echo "Add package $request"
                echo "ADD_SUBDIRECTORY(packages/${request})" >> CMakeLists.txt
                One=True
                Any=True
                break
            fi
        done
        if [ $One == "False" ] ; then
            echo "${RED}WARNING: ${request} is not valid package. Make sure it is in packages/*!${NOCOLOR}"
        fi
    done

    if [ $Any == "False" ] ; then
        echo "${RED}WARNING: No valid package is specified. Make sure the names are consistent with packages/*!${NOCOLOR}"
        exit
    else
        cp CMakeLists.txt ..
    fi
fi

# For different servers, load different modules
if [ "$machine" == "Linux" ] 
then
    case "$HOSTNAME" in
        sh-ln*.stanford.edu) 
            module load cmake
            module load math
            module load py-tensorflow/1.9.0_py36
            module load py-keras/2.2.4_py36
            module load cuda/9.0.176
            module load system
            module load x11
            module load openblas
            ;;
        *)  echo "${RED}Unknown linux machine!${NOCOLOR}" ;;
    esac
fi
    
if [ "$machine" == "Mac" ] ;
then
    cmake -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++)  -DPYTHON_INCLUDE_DIR=$(python3 -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")  -DPYTHON_LIBRARY=$(python3 -c "import distutils.sysconfig as sysconfig; print(sysconfig.get_config_var('LIBDIR'))") ..
fi

if [ "$machine" == "Linux" ] ;
then
    cmake -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++) -DPYTHON_INCLUDE_DIR=$(python3 -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")  -DPYTHON_LIBRARY=$(python3 -c "import distutils.sysconfig as sysconfig; print(sysconfig.get_config_var('LIBDIR'))")  ..
fi

make

cd ..
#rm -rf build
rm CMakeLists.txt
