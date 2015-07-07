[![Build Status](https://travis-ci.org/jbakosi/quinoa.svg?branch=master)](https://travis-ci.org/jbakosi/quinoa)
[![Issues](https://img.shields.io/github/issues/jbakosi/quinoa.svg)](https://github.com/jbakosi/quinoa/issues)
[![License](https://img.shields.io/github/license/jbakosi/quinoa.svg)](#license)
<!---[![Release](https://img.shields.io/github/release/jbakosi/quinoa.svg)](https://github.com/jbakosi/quinoa/releases/latest)-->

## What is Quinoa?

Quinoa is a set of computational tools that enables research and numerical
analysis in fluid dynamics. At this time it is a test-bed to experiment with
various algorithms using fully asynchronous runtime systems.

## Organization

Currently, Quinoa consists of the following tools:
  - [<B>inciter</B>](https://jbakosi.github.io/quinoa/inciter_doc.html) - Solve
    a PDE on an unstructured mesh using overdecomposition
  - [<B>walker</B>](https://jbakosi.github.io/quinoa/walker_doc.html) - Random
    walker for stochastic differential equations
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

See the [LICENSE](https://github.com/jbakosi/quinoa/blob/master/LICENSE).

## Authors

Jozsef Bakosi (jbakosi@lanl.gov)

## How to build

### 1. Pick compilers

   ```
    $ CC=mpicc CXX=mpic++ FC=mpif90
   ```
   - Currently, MPI is required, use the OpenMPI wrappers
   - The underlying C++ compiler must support the C++11 standard

### 2. Build the third-party libraries

   ```
    $ cd <quinoa>/tpl; mkdir build; cd build
    $ cmake ..
    $ make
   ```

### 3. Build Quinoa

   ```
    $ cd <quinoa>; mkdir build; cd build
    $ cmake ../src
    $ make
   ```
   - All executables will be in <tt>./Main</tt>
