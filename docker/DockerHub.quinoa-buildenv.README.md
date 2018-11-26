## Quinoa

_Adaptive computational fluid dynamics_ - https://quinoacomputing.github.io

Quinoa is a set of computational tools that enables research and numerical analysis in fluid dynamics. Using the [Charm++](http://charmplusplus.org) runtime system, we employ _asynchronous_ (or non-blocking) parallel programming and decompose computational problems into a large number of work units (that may be more than the available number of processors) enabling _arbitrary overlap_ of parallel computation, communication, input, and output. Then the runtime system _dynamically_ and _automatically_ homogenizes computational load across the simulation distributed across many computers.

Our ultimate goal is to simulate large and complex engineering multiphysics problems with a production-quality code that is extensible and maintainable, using hardware resources efficiently, even for problems with _a priori_ unknown, heterogeneous, and dynamic load distribution.

For more details on philosophy, documentation, software design, journal papers, license, contributing see the
[documentation](https://quinoacomputing.github.io).

### How to use these images
The images in this repository are configured for a single computer, and thus intended to quickly try the pre-built executables on a multi-core workstation. For production runs on clusters of networked compute nodes you should build from source. See https://quinoacomputing.github.io for more details.

Run the unit-, and regression tests in the Alpine container on your local machine:

    docker run -ti quinoacomputing/quinoa:alpine
    cd quinoa/build && ../script/run_tests.sh
