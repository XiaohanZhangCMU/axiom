# axiom: gpu/c++ & python binding playground using boost and pybind11.
Add C++ program and simulation packages to axiom/ and export the variables and methods to python. Ideally, many simulation models can be glued together in this way.  

System tested:
1) OS X > Yosemite 2) Linux: CentOS >= 7.5.1804 or ubuntu >=16.04.4 LTS).
3) python: 2.7 or >=3. 
4) boost = 1.66.0, 1.68.0.

Install:

0) git clone git@github.com:XiaohanZhangCMU/axiom.git
1) cd axiom; 
2) choose desired examples (see below) to build in CMakeLists.txt::EXAMPLES part, for example, ADD_SUBDIRECTORY(mdsw)  
3) mkdir build; cd build; 
4) cmake ..
5) make

Optional packages: (required for mdsw and mdfem)
6) pip install PyOpenGL PyOpenGL_accelerate 
7) python3 -m pip install Pillow
Main functionals of axiom:

Examples:
1) zoo: playground of simple pybind11 binding. 
2) tensor: gpu and boost::share_ptr for automatic gpu memeory allocation and syncronization. 
2) mdsw and mdfem. two modules from MD++ (http://micro.stanford.edu/MDpp), a molecular dynamics simulation package written in C++. 

 

Bug reports or questions:

xiaohanzhang.cmu@gmail.com
caiwei@stanford.edu


