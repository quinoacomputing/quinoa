## Quinoa

<img src="https://quinoacomputing.github.io/quinoa.svg" align="right" width="25%" background=transparent>

_Adaptive computational fluid dynamics_ - https://quinoacomputing.github.io

Quinoa is a set of computational tools that enables research and numerical
analysis in fluid dynamics. Using the [Charm++](http://charmplusplus.org)
runtime system, we employ _asynchronous_ (or non-blocking) parallel programming
and decompose computational problems into a large number of work units (that may
be more than the available number of processors) enabling _arbitrary
overlap_ of parallel computation, communication, input, and output. Then the
runtime system _dynamically_ and _automatically_ homogenizes computational load
across the simulation distributed across many computers.

Our ultimate goal is to simulate large and complex engineering multiphysics
problems with a production-quality code that is extensible and maintainable,
using hardware resources efficiently, even for problems with _a priori_
unknown, heterogeneous, and dynamic load distribution.

## Directory layout

    quinoa
      ├── cmake/ - CMake code shared between github.com/quinoacomputing/quinoa-tpl and src/.
      ├── doc/ - Documentation, rendered at quinoacomputing.org.
      ├── external/ - External packages pulled in as submodules from github.com/quinoacomputing/quinoa-tpl.
      ├── src/ - Compilable sources. For details see quinoacomputing.org/files.html.
      ├── tests/ - Unit-, and regression tests.
      ├── tools/ - Development utilities and docker files.
      ├── LICENSE - Copyright and license.
      └── README.md - This file, rendered at github.com/quinoacomputing/quinoa.

## More info

For more details on philosophy, documentation, software design, journal papers,
license, and contributing see the
[documentation](https://quinoacomputing.github.io).

