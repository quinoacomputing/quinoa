## What is Quinoa?

<img src="https://quinoacomputing.github.io/quinoa/quinoa.svg" align="right" width="25%" background=transparent>
Quinoa is a set of computational tools that enables research and numerical analysis in fluid dynamics. At this time it is a test-bed to experiment with various algorithms using fully asynchronous runtime systems.

## Organization

Currently, Quinoa consists of the following tools:
  - [<B>inciter</B>](http://quinoacomputing.github.io/quinoa/inciter_doc.html) - Solve a system of partial differential equations discretized on a 3D unstructured mesh using overdecomposition
  - [<B>walker</B>](http://quinoacomputing.github.io/quinoa/walker_doc.html) - Random walker for stochastic differential equations
  - [<B>rngtest</B>](http://quinoacomputing.github.io/quinoa/rngtest_doc.html) - Random number generator test suite
  - [<B>unittest</B>](http://quinoacomputing.github.io/quinoa/unittest_doc.html) - Unit test suite
  - [<B>meshconv</B>](http://quinoacomputing.github.io/quinoa/meshconv_doc.html) - Mesh file converter

## Build

### 1. Install prerequisites

- Debian/Ubuntu linux: (line 1: required, line 2: recommended)

   ```
   apt-get install cmake gfortran gcc g++ openmpi-bin libopenmpi-dev
   apt-get install gmsh libpugixml-dev libpstreams-dev libboost-all-dev liblapack-dev liblapacke-dev libhdf5-dev libhdf5-openmpi-dev libhypre-dev
   ```

- Mac OS X: (line 1: required, line 2: recommended)

   ```
   port install cmake openmpi-clang38 && port select clang mp-clang-3.8 && port select mpi openmpi-clang38-fortran
   port install gmsh pugixml boost hdf5 +hl +openmpi hypre +openmpi
   ```

### 2. Clone, build third-party libraries, build & test

   ```
   git clone https://github.com/quinoacomputing/quinoa.git; cd quinoa
   mkdir tpl/build; cd tpl/build; cmake ..; make; cd -
   mkdir build; cd build; cmake ../src; make; ../script/run_tests.sh
   ```

   - All executables will be in <tt>quinoa/build/Main</tt>

## Authors

Jozsef Bakosi (jbakosi@lanl.gov)

## Resources

 - [<B>Source</B>](https://github.com/quinoacomputing/quinoa)
 - [<B>License</B>](https://github.com/quinoacomputing/quinoa/blob/master/LICENSE)
 - [<B>Documentation</B>](http://quinoacomputing.github.io/quinoa/index.html)
 - [<B>Manifesto</B>](http://quinoacomputing.github.io/quinoa/why.html)
 - [<B>Roadmap</B>](https://github.com/quinoacomputing/quinoa/issues)
