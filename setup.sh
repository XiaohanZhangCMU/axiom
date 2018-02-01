#!/bin/sh

if [ ! -d build ]; then
    mkdir build
fi
cd build

# Check which system you are running

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${unameOut}"
esac

# use gcc
# CC=$(command -v gcc)
# CXX=$(command -v g++)

# use icc
CC=$(command -v icc)
CXX=$(command -v icpc)

if [ ${machine} = "Mac" ]; then 
    if [ ! -d ${HOME}/boost ]; then
        echo "Boost has been installed to ${HOME}/boost. Skip install!"
    else
        echo "Boost is missing from ${HOME}/boost. Installing boost. Patience!!!"

        if type xcode-select >&- && xpath=$( xcode-select --print-path ) && test -d "${xpath}" && test -x "${xpath}" ; then
            echo "xcode-select is correctly installed"
        else
            echo "xcode-select isn't correctly installed. Install now!"
            xcode-select --install
        fi

        if ! [ -x "$(command -v brew)" ]; then
            echo 'Error: brew needs to be installed.' >&2
            xcode-select --install
            ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
            brew doctor
        fi
        brew install wget --with-libressl
        brew upgrade gcc
        pushd /usr/local/bin
        ln -s g++-7 g++
        ln -s gcc-7 gcc
        pushd
        export PATH="/usr/local/bin:\$PATH"
        echo "export PATH=\"/usr/local/bin:\$PATH\"" >> ${HOME}/.bash_profile

        brew install openblas

        miniconda3_install_sh="Miniconda3-latest-MacOSX-x86_64.sh"

        if ! [ -x "$(command -v conda)" ]; then
            echo 'Error: miniconda3 needs to be installed with root as ${HOME}/miniconda3' >&2
            if [ ! -f $miniconda3_install_sh ]; then
                wget "https://repo.continuum.io/miniconda/${miniconda3_install_sh}"
            fi
            chmod +x ${miniconda3_install_sh}
            ./${miniconda3_install_sh}
            echo "Switch back python. Py > 3.6.0 leads to strange runtime error for import_array()"
            conda install python=3.6.0 
            conda install numpy
        fi
        export CPLUS_INCLUDE_PATH="${HOME}/miniconda3/include/python3.6m:${CPLUS_INCLUDE_PATH}"
        export CPLUS_INCLUDE_PATH="${HOME}/miniconda3/lib/python3.6/site-packages/numpy/core/include:${CPLUS_INCLUDE_PATH}"
        echo "export CPLUS_INCLUDE_PATH=\"\${HOME}/miniconda3/include/python3.6m:\${CPLUS_INCLUDE_PATH}\"" >> ${HOME}/.bash_profile
        echo "export CPLUS_INCLUDE_PATH=\"\${HOME}/miniconda3/lib/python3.6/site-packages/numpy/core/include:\${CPLUS_INCLUDE_PATH}\"" >> ${HOME}/.bash_profile

        boost_install="boost_1_66_0"
        if [ ! -f "${boost_install}.tar.gz" ]; then
            wget "https://dl.bintray.com/boostorg/release/1.66.0/source/${boost_install}.tar.gz"
        fi
        if [ ! -d ${boost_install} ]; then
            tar -xf ${boost_install}.tar.gz
        fi
        cd ${boost_install}
        ./bootstrap.sh -with-libraries=python -with-libraries=thread --prefix=${HOME}/boost/
        ./b2 install 
        export DYLD_LIBRARY_PATH=${HOME}/boost/lib/:${DYLD_LIBRARY_PATH}
    fi
fi


cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX ..
make
#python ../test_axiom.py
cd ..

