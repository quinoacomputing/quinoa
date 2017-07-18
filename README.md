## What is Quinoa?

<img src="https://quinoacomputing.github.io/quinoa.svg" align="right" width="25%" background=transparent>
Quinoa is a set of computational tools that enables research and numerical analysis in fluid dynamics. At this time it is a test-bed to experiment with various algorithms using fully asynchronous runtime systems.

## Organization

Currently, Quinoa consists of the following tools:
  - [<B>inciter</B>](http://quinoacomputing.github.io/inciter_doc.html) - Partial differential equations solver on 3D unstructured grids
  - [<B>walker</B>](http://quinoacomputing.github.io/walker_doc.html) - Random walker for stochastic differential equations
  - [<B>rngtest</B>](http://quinoacomputing.github.io/rngtest_doc.html) - Random number generator test suite
  - [<B>unittest</B>](http://quinoacomputing.github.io/unittest_doc.html) - Unit test suite
  - [<B>meshconv</B>](http://quinoacomputing.github.io/meshconv_doc.html) - Mesh file converter

## Build

#### 1. Install prerequisites

- Debian/Ubuntu linux: (line 1: required, line 2: recommended)

   ```
   apt-get install cmake gfortran gcc g++ openmpi-bin libopenmpi-dev
   apt-get install gmsh libpugixml-dev libpstreams-dev libboost-all-dev liblapack-dev liblapacke-dev libhdf5-dev libhdf5-openmpi-dev libhypre-dev
   ```

- Mac OS X: (line 1: required, line 2: recommended)

   ```
   port install cmake openmpi-clang38 && port select clang mp-clang-3.8 && port select mpi openmpi-clang38-fortran
   port install gmsh pugixml lapack boost
   ```

#### 2. Clone, build third-party libraries, build & test

   ```
   git clone --recursive https://github.com/quinoacomputing/quinoa.git; cd quinoa
   mkdir tpl/build; cd tpl/build; cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 ..; make; cd -
   mkdir build; cd build; cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc ../src; make; ../script/run_tests.sh
   ```

   - All executables will be in <tt>quinoa/build/Main</tt>

## Authors
(in chronological order of contribution)

1. [Jozsef Bakosi](https://github.com/jbakosi)
2. [Francisco Gonzalez](https://github.com/franjgonzalez)
3. [Brandon Rogers](https://github.com/brog2610)
4. [Christoph Junghans](https://github.com/junghans)
5. [Robert Pavel](https://github.com/rspavel)
5. [Aditya Pakki](https://github.com/adityapakki)

## Resources

[source](https://github.com/quinoacomputing/quinoa) - [license](https://github.com/quinoacomputing/quinoa/blob/master/LICENSE) - [documentation](http://quinoacomputing.github.io/index.html) - [manifesto](http://quinoacomputing.github.io/why.html) - [roadmap](https://github.com/quinoacomputing/quinoa/issues) - [docker](https://hub.docker.com/r/quinoacomputing) - [travis](https://travis-ci.org/quinoacomputing/quinoa) - [codecov](https://codecov.io/gh/quinoacomputing/quinoa/commits) - [gcov](http://quinoacomputing.github.io/coverage.html) - [cppcheck](http://quinoacomputing.github.io/cppcheck/index.html) - [doxygen](http://quinoacomputing.github.io) - [sonarqube](https://sonarqube.com/organizations/quinoacomputing) - [analytics](https://www.openhub.net/p/quinoacomputing) - [cite](https://zenodo.org/badge/latestdoi/38454430)
