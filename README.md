# evo_q
A simple library for using genetic algorithms with c++ and python

# Installation
To install, clone the repository into a known location. Change directory to the installed path and run the following commands:

    mkdir build && cd build
    cmake ..
    make

This will generate the required shared libraries in the build directory. These libraries can then be copied into your project and linked.

To install the python bindings run

    git clone https://github.com/sam-stromswold/evo_q
    pip install ./evo_q
    
You may then use the library in your packages by including 

    import evo_p
    
In the beginning of your python scripts.
