[![Join the chat at https://gitter.im/jbakosi/quinoa](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/jbakosi/quinoa?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## What is Quinoa?

Quinoa is a set of computational tools that enables research and numerical analysis in fluid dynamics. At this time it is a test-bed to experiment with various algorithms using fully asynchronous runtime systems.

## Organization

Currently, Quinoa consists of the following tools:
  - [<B>inciter</B>](https://jbakosi.github.io/quinoa/inciter_doc.html) - Solve a PDE on an unstructured mesh using overdecomposition
  - [<B>walker</B>](https://jbakosi.github.io/quinoa/walker_doc.html) - Random walker for stochastic differential equations
  - [<B>rngtest</B>](https://jbakosi.github.io/quinoa/rngtest_doc.html) - Random number generator test suite
  - [<B>unittest</B>](https://jbakosi.github.io/quinoa/unittest_doc.html) - Unit test suite
  - [<B>meshconv</B>](https://jbakosi.github.io/quinoa/meshconv_doc.html) - Mesh file converter

## Goals

  - [<B>Designed for the exascale era</B>](https://jbakosi.github.io/quinoa/why.html#exascale)
  - [<B>Verified and proven to be correct</B>](https://jbakosi.github.io/quinoa/why.html#correct)
  - [<B>Optimized for performance, power, and reliability</B>](https://jbakosi.github.io/quinoa/why.html#optimized)
  - [<B>Advanced computer science outsourced to experts</B>](https://jbakosi.github.io/quinoa/why.html#outsource)
  - [<B>Using a language that can cope with complexity</B>](https://jbakosi.github.io/quinoa/why.html#language)
  - [<B>Highly-valued programmer productivity</B>](https://jbakosi.github.io/quinoa/why.html#productivity)
  - [<B>User and developer friendly</B>](https://jbakosi.github.io/quinoa/why.html#friendly)
  - [<B>Well documented</B>](https://jbakosi.github.io/quinoa/why.html#documented)
  - [<B>Fun to work on</B>](https://jbakosi.github.io/quinoa/why.html#fun)

## License

See the [license](https://github.com/jbakosi/quinoa/blob/master/LICENSE)

## Authors

Jozsef Bakosi (jbakosi@lanl.gov)

## Build on linux or mac

### 1. Install prerequisites

#### Debian/Ubuntu:
   ```
    $ apt-get install git cmake gfortran gcc g++ gmsh libpugixml-dev libpstreams-dev libboost-all-dev liblapack-dev liblapacke-dev libhdf5-dev libhdf5-openmpi-dev libhypre-dev
   ```
#### Mac OS X:
   ```
    $ port install git cmake openmpi-clang38 gmsh pugixml boost
    $ port install hdf5 +hl +openmpi
    $ port install hypre +openmpi
   ```

### 2. Clone

   ```
    $ git clone https://github.com/jbakosi/quinoa.git
   ```

### 3. Build the third-party libraries

   ```
    $ CC=mpicc CXX=mpic++ FC=mpif90
    $ cd quinoa/tpl; mkdir build; cd build
    $ cmake ..
    $ make -sj4
   ```

### 4. Build & test Quinoa

   ```
    $ cd quinoa; mkdir build; cd build
    $ cmake ../src
    $ make -sj4
    $ Main/charmrun +p4 Main/unittest
    $ ctest -j4 -LE stringent
   ```
   - All executables will be in <tt>quinoa/build/Main</tt>
   - Browse the [documentation](http://jbakosi.github.io/quinoa/index.html)
