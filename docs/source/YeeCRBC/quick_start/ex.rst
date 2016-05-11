Running the built-in Yee shceme CRBC examples
=============================================

Install cmake (version >=3.0), if you don't have one. ::

  cmake -version

By default, 'sh cmake-3.3.2-Linux-x86_64.sh' can install at current dir. Then create a link to cmake and put in the ~/bin and add this to PATH variable.

Setup variable CC, CXX to indicate which compiler one want to use to compile the CRBC yee library(i.e. gcc or icc, g++ or icpc, or by default). Note that this have to be consistent if one want to manually link the library to their own yee scheme code later. In csh shell,::

  setenv CC gcc
  setenv CXX g++

By default yee CRBC library will be install to */usr/local/lib* and header files in */usr/local/include/yeecrbc*. If you don't have access to those or want to install it locally:: 

  tar -zxvf yee_crbc_lib.tar.gz
  mkdir build_gcc
  cd build_gcc
  cmake ../yee_crbc_lib_1.0.0 -DCMAKE_INSTALL_PREFIX:PATH=~/YeeCRBC
  make
  make install

For some compiler, one may need to add flag "-DCMAKE_CXX_FLAGS:STRING=-std=c++11" to the cmake command. In this case, the library and header files will be placed at *~/YeeCRBC/lib* and *~/YeeCRBC/include/crbc* repectively. To remove the build, just delete those directories::

  cd ~/YeeCRBC
  rm -rf build_gcc lib include

To test the build, go to 2d_yee(4500 iterations) or 3d_yee(1200 iterations), e.g.::

  cd build_gcc/examples/2d_yee
  ./yee_TM.x

To build the example again(e.g. 2d_yee),::

  cd build_gcc/examples/2d_yee
  setenv yeecrbc_DIR ~/YeeCRBC/lib
  cmake .
  make

To manually compile and link the example(e.g. 2d_yee), note that one need to specify the environment variable $yeecrbc_DIR::

  gcc -L$yeecrbc_DIR -I$yeecrbc_DIR/../include/crbc yee_TM.c -o yee_TM.x -lyeecrbc -lstdc++ -lm

or::

  g++ -L$yeecrbc_DIR -I$yeecrbc_DIR/../include/crbc yee_TM.c -o yee_TM.x -lyeecrbc

See User Guides for more compiling options and running examples.

