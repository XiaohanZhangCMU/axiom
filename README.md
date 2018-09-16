# openAI/gym environment for Molecular Dynamics 

Why called axiom? **aXioM** == **Xiao** Han Zhang + **Mo** Han Jiang :)

### Axiom is a framework that allows you to quickly bind a GPU-accelerated C++ project with python-3 using pybind11 (https://github.com/pybind/pybind11, a wrapper of boost)  

### For scientific computing comunity, it is ideal that one can move the computation laborous part of your code to C++, and call them from python to interface many other available packages. A direct benefit of using axiom is making your own C++ simulator an OpenAI/gym environment and run all kinds of Reinforcement learning algorithms easily.  


Steps to follow to create your own environment:
1) Add your GPU-C++ project to axiom/ as a new subdirectory. 
2) Write your own CMakeLists.txt in the subdirectory. Also, write your own _axiom.cpp in a similar way as _axiom.cpp in the other examples subdirectories. This exports the variables and methods you desired to python  
3) CMakeLists.txt of the root directory compiled all subdirectories in the EXAMPLE section to generate a dynamic library with the same name as the subdirectory and put in axiom/lib.   
4) In python script, you need to add sys.path.append('your_path_to/axiom/lib/') at the top of the python script. 
5) Program the interface functions required by gym (https://github.com/openai/gym/wiki/Environments).
6) Enjoy!

## Features:
1) Implement greedy search algorithm for Frank Partial Dislocation energy barrier calculation in silicon thin film.
2) On going: Deep Q network search algorithm for (0).
3) tensor module has realized automatic memory synchronization between GPU and CPU for tensors.
4) A view helper script is provided as scripts/tests/View.py using PyOpenG. You can develop similar one for your own simulator. .
5) Export C++ simulation packages to openAI/gym using this framework. Examples: mdsw and mdfem.
6) Arrays can be exported as numpy arrays and modified in place (no deep memory copy is needed).

## System and software tested:
1) OS X > Yosemite, Linux: CentOS >= 7.5.1804 or ubuntu >=16.04.4 LTS).
2) python: 2.7 or >=3. 
3) boost = 1.66.0, 1.68.0.

## Install mdsw and compile:

1) git clone git@github.com:XiaohanZhangCMU/axiom.git; cd axiom;   
   in the CMakeLists.txt::EXAMPLES part, add ADD_SUBDIRECTORY(mdsw). 
2) mkdir build; cd build; 
3) cmake ..
4) make  
Optional packages for mdsw and mdfem. Required only if you want to visualize mdsw.   
5) pip install PyOpenGL PyOpenGL_accelerate (On Ubuntu or Linux, you need to do: sudo apt-get install python-opengl as well)
6) python3 -m pip install Pillow   

## Examples:
0) zoo: playground of simple pybind11 binding. 
1) tensor: gpu and boost::share_ptr for automatic gpu memeory allocation and syncronization. 
2) mdsw and mdfem. two modules from MD++ (http://micro.stanford.edu/MDpp), a molecular dynamics simulation package written in C++. 

<img src="animation.gif" alt="Drawing" style="width: 400px;"/>

Bug reports or questions:

xiaohanzhang.cmu@gmail.com


