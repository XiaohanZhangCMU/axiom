# axiom: gpu/c++ & python binding playground using boost and pybind11.
0) axiom uses pybind11 (https://github.com/pybind/pybind11) a wrapper of boost to bind. 
1) Add C++/GPU program and simulation packages to axiom/ as a new subdirectory and export the variables and methods to python. This requires you write your own CMakeLists.txt in the subdirectory. Also, write your own _axiom.cpp in a similar way as _axiom.cpp in the other examples subdirectories.
2) Ideally, many simulation models can be glued together in this framework.  
3) CMakeLists.txt of the root directory compiled all subdirectories in the EXAMPLE section to generate a dynamic library with the same name as the subdirectory and put in axiom/lib.   
4) In python script, you need to add sys.path.append('your_path_to/axiom/lib/') at the top of the python script. 
5) mdsw and mdfem are two modules of MD++ (http://micro.stanford.edu/MDpp) and are combined as a part of openAI:gym. 
6) View.py is a visulization helper script. You can develop similar one for your own simulator. 


## System tested:
1) OS X > Yosemite 2) Linux: CentOS >= 7.5.1804 or ubuntu >=16.04.4 LTS).
3) python: 2.7 or >=3. 
4) boost = 1.66.0, 1.68.0.

## Install:

0) git clone git@github.com:XiaohanZhangCMU/axiom.git
1) cd axiom; 
2) choose desired examples (see below) to build in CMakeLists.txt::EXAMPLES part, for example, ADD_SUBDIRECTORY(mdsw)  
3) mkdir build; cd build; 
4) cmake ..
5) make

Optional packages for mdsw and mdfem. Required if you want to visualize atoms. 
6) pip install PyOpenGL PyOpenGL_accelerate 
7) python3 -m pip install Pillow

## Examples:
1) zoo: playground of simple pybind11 binding. 
2) tensor: gpu and boost::share_ptr for automatic gpu memeory allocation and syncronization. 
2) mdsw and mdfem. two modules from MD++ (http://micro.stanford.edu/MDpp), a molecular dynamics simulation package written in C++. 

 
Bug reports or questions:

xiaohanzhang.cmu@gmail.com
caiwei@stanford.edu


