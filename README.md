<!--[![Travis Build Status](https://travis-ci.org/jbakosi/quinoa.svg?branch=master)](https://travis-ci.org/jbakosi/quinoa)-->
[![Drone Build Status](http://bakosi.com:8000/api/badge/github.com/jbakosi/quinoa/status.svg?branch=master)](http://bakosi.com:8000/github.com/jbakosi/quinoa)
[![Issues](https://img.shields.io/github/issues/jbakosi/quinoa.svg)](https://github.com/jbakosi/quinoa/issues)
[![Join the chat at https://gitter.im/jbakosi/quinoa](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/jbakosi/quinoa?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![License](https://img.shields.io/github/license/jbakosi/quinoa.svg)](https://github.com/jbakosi/quinoa/blob/master/LICENSE)

<!---[![Release](https://img.shields.io/github/release/jbakosi/quinoa.svg)](https://github.com/jbakosi/quinoa/releases/latest)-->

## What is Quinoa?

Quinoa is a set of computational tools that enables research and numerical analysis in fluid dynamics. At this time it is a test-bed to experiment with various algorithms using fully asynchronous runtime systems.

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

## Authors

Jozsef Bakosi (jbakosi@lanl.gov)

## Get started

### 1. Install prerequisites
   Click on the links below for details on how to install the prerequisites for
   your system.
   - [debian/testing/gnu](https://github.com/jbakosi/quinoa-docker/blob/master/debian/testing/gnu/alltpl/Dockerfile)
   - [debian/testing/clang](https://github.com/jbakosi/quinoa-docker/blob/master/debian/testing/clang/alltpl/Dockerfile)
   - [ubuntu/precise/gnu](https://github.com/jbakosi/quinoa-docker/blob/master/ubuntu/precise/gnu/alltpl/Dockerfile)
   - [ubuntu/trusty/gnu](https://github.com/jbakosi/quinoa-docker/blob/master/ubuntu/trusty/gnu/alltpl/Dockerfile)
   - [add yours](https://github.com/jbakosi/quinoa/issues/73)

### 2. Clone

   ```
    $ git clone https://github.com/jbakosi/quinoa.git
   ```

### 3. Pick compilers

   ```
    $ CC=mpicc CXX=mpic++ FC=mpif90
   ```
   - Currently, MPI is required, use the OpenMPI wrappers
   - The underlying C++ compiler must support the C++11 standard

### 4. Build the third-party libraries

   ```
    $ cd quinoa/tpl; mkdir build; cd build
    $ cmake ..
    $ make
   ```

### 5. Build Quinoa

   ```
    $ cd quinoa; mkdir build; cd build
    $ cmake ../src
    $ make
   ```
   - All executables will be in <tt>quinoa/build/Main</tt>

### 6. Contribute
   - Browse the [documentation](http://jbakosi.github.io/quinoa/index.html)
   - Check out the issues labeled [help wanted](https://github.com/jbakosi/quinoa/labels/help%20wanted)
