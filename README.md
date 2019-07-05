# evo_q
A simple library for using genetic algorithms with c++ and python

# Installation
To install, clone the repository into a known location. Change directory to the installed path and run the following commands:

    mkdir build && cd build
    cmake ..
    make

This will generate the required shared libraries in the build directory. These libraries can then be copied into your project and linked. To use the python bindings copy the shared library starting with evo_p into your working directory and link using

    import evo_p
