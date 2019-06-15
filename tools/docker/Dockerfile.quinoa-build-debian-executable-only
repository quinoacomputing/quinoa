################################################################################
# vim: filetype=dockerfile:
#
# \file      tools/docker/Dockerfile.quinoa-build-debian-executabl-only
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Dockerfile for a testing a single executable of Quinoa on Debian
# \details   Dockerfile for a building and testing a single executable of
#            Quinoa on Debian Linux using the native system-wide compiler and
#            MPI wrappers.
#
# We start from a Debian Linux distribution and install all possible
# prerequisites, then build the required but missing third-party libraries,
# build the code, and test.
#
# This docker image is intended to test building a single executable, preceded
# by cloning and building only those third-party libraries that are required
# for the executable.
#
################################################################################

FROM debian:stable

# Install system-wide prerequisites
RUN apt-get update -y && apt-get install -y m4 autoconf git cmake gfortran gcc g++ openmpi-bin libopenmpi-dev gmsh libpugixml-dev libpstreams-dev libboost-all-dev libblas-dev liblapack-dev liblapacke-dev zlib1g-dev libhdf5-dev libhdf5-openmpi-dev binutils-dev libx11-dev libxpm-dev libxft-dev libxext-dev ninja-build flex bison libdw-dev libdwarf-dev vim wget

# Setup user
RUN adduser --gecos "" --disabled-password quinoa
USER quinoa
WORKDIR /home/quinoa
CMD ["/bin/bash"]
SHELL ["/bin/bash", "-c"]

# Clone quinoa non-recursively
RUN git clone http://github.com/quinoacomputing/quinoa.git
# Checkout commit to be tested
ARG COMMIT
RUN cd quinoa && git checkout $COMMIT && git log -1 HEAD
# Update submodules required by executable tested
RUN wget https://raw.githubusercontent.com/quinoacomputing/quinoa/$COMMIT/doc/pages/build.dox
ARG EXECUTABLE
ENV EXECUTABLE $EXECUTABLE
# Get git update submodule command for executable
RUN grep "^@ref $EXECUTABLE" build.dox -A2 | grep -v $EXECUTABLE | grep -v -e "^$" > submodule_update.sh
RUN cat submodule_update.sh
# Pull in submodules required for the executable
RUN cd quinoa && git submodule init && git submodule update && cd external && sh ../../submodule_update.sh && cd .. && git submodule status --recursive
# Build TPLs
RUN cd quinoa && mkdir -p external/build && cd external/build && cmake -D${EXECUTABLE^^}_ONLY=true -DCMAKE_BUILD_TYPE=Release .. && make -sj$(grep -c processor /proc/cpuinfo)
# Build code
RUN cd quinoa && mkdir -p build && cd build && cmake -DCMAKE_CXX_FLAGS=-Werror -DCMAKE_BUILD_TYPE=Release -GNinja -DRUNNER=mpirun -DRUNNER_NCPUS_ARG=-n -DRUNNER_ARGS=-oversubscribe ../src && ninja ${EXECUTABLE,,}
# Run tests
RUN cd quinoa/build && npe=$(grep -c processor /proc/cpuinfo) && if [ ${EXECUTABLE} = unittest ]; then mpirun -n $npe Main/unittest -v -q; else ctest -j $npe --output-on-failure -LE extreme -R ${EXECUTABLE,,}; fi
