Quinoa
======

_Adaptive computational fluid dynamics_

Quinoa is a set of computational tools that enables research and numerical analysis in fluid dynamics. At this time it is a test-bed to experiment with various algorithms using fully asynchronous runtime systems.

Organization
============

Quinoa consists of the following tools:

  - [inciter](http://quinoacomputing.github.io/inciter_doc.html) - Partial differential equations solver on 3D unstructured grids
  - [walker](http://quinoacomputing.github.io/walker_doc.html) - Random walker for stochastic differential equations
  - [rngtest](http://quinoacomputing.github.io/rngtest_doc.html) - Random number generator test suite
  - [unittest](http://quinoacomputing.github.io/unittest_doc.html) - Unit test suite
  - [meshconv](http://quinoacomputing.github.io/meshconv_doc.html) - Mesh file converter

Build
=====

1. Install prerequisites

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

2. Clone, build third-party libraries, build & test

```
git clone --recursive https://github.com/quinoacomputing/quinoa.git; cd quinoa
mkdir tpl/build; cd tpl/build; cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 ..; make; cd -
mkdir build; cd build; cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc ../src; make; ../script/run_tests.sh
```

All executables will be in <tt>quinoa/build/Main</tt>

Authors
=======
(in chronological order of contribution)

1. [Jozsef Bakosi](https://github.com/jbakosi)
2. [Francisco Gonzalez](https://github.com/franjgonzalez)
3. [Brandon Rogers](https://github.com/brog2610)
4. [Christoph Junghans](https://github.com/junghans)
5. [Robert Pavel](https://github.com/rspavel)
6. [Aditya Pakki](https://github.com/adityapakki)
7. [Bob Bird](https://github.com/rfbird)
8. [Aditya Pandare](https://github.com/adityakpandare)

Resources
=========

- Source code --- https://github.com/quinoacomputing/quinoa
- License --- https://github.com/quinoacomputing/quinoa/blob/master/LICENSE
- Documentation --- http://quinoacomputing.github.io/index.html
- Manifesto --- http://quinoacomputing.github.io/why.html
- Roadmap --- https://github.com/quinoacomputing/quinoa/issues
- Docker --- https://hub.docker.com/r/quinoacomputing
- Travis --- https://travis-ci.org/quinoacomputing/quinoa
- Codecov --- https://codecov.io/gh/quinoacomputing/quinoa/commits
- Gcov --- http://quinoacomputing.github.io/coverage.html
- CppCheck --- http://quinoacomputing.github.io/cppcheck/index.html
- Doxygen --- http://quinoacomputing.github.io
- Sonarqube --- https://sonarqube.com/organizations/quinoacomputing
- Analytics --- https://www.openhub.net/p/quinoacomputing
- Cite --- https://zenodo.org/badge/latestdoi/38454430
- CLA --- https://www.clahub.com/agreements/quinoacomputing/quinoa

License
=======

> Copyright (c) 2016-2018, Los Alamos National Security, LLC
> All rights reserved.
> 
> Copyright 2016-2018. Los Alamos National Security, LLC. This software was
> produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos
> National Laboratory (LANL), which is operated by Los Alamos National Security,
> LLC for the U.S. Department of Energy. The U.S. Government has rights to use,
> reproduce, and distribute this software. NEITHER THE GOVERNMENT NOR LOS ALAMOS
> NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
> LIABILITY FOR THE USE OF THIS SOFTWARE. If software is modified to produce
> derivative works, such modified software should be clearly marked, so as not to
> confuse it with the version available from LANL.
>  
> Additionally, redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are met:
> 
> 1. Redistributions of source code must retain the above copyright notice, this
> list of conditions and the following disclaimer.
> 
> 2. Redistributions in binary form must reproduce the above copyright notice,
> this list of conditions and the following disclaimer in the documentation and/or
> other materials provided with the distribution.
> 
> 3. Neither the name of Los Alamos National Security, LLC, Los Alamos National
> Laboratory, LANL, the U.S. Government, nor the names of its contributors may be
> used to endorse or promote products derived from this software without specific
> prior written permission.
> 
> THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
> "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
> THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
> ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR
> CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
> OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
> SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
> INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
> CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
> IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
> OF SUCH DAMAGE.
