################################################################################
#
# \file      script/build_openmpi.sh
# \author    J. Bakosi
# \date      Tue 07 Jul 2015 06:19:16 AM MDT
# \copyright 2012-2015, Jozsef Bakosi.
# \brief     Download and build OpenMPI
#
################################################################################

#!/bin/bash

sudo apt-get install -y clang-3.4
wget http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.6.tar.gz
tar xzf openmpi-1.8.6.tar.gz
cd openmpi-1.8.6
OPENMPI=$HOME/openmpi
./configure --prefix=$OPENMPI FC=gfortran
make -sj2 install
cd -

export PATH=$OPENMPI/bin:$PATH
export LD_LIBRARY_PATH=$OPENMPI/lib:$LD_LIBRARY_PATH
export CC=$OPENMPI/bin/mpicc
export CXX=$OPENMPI/bin/mpic++
export FC=$OPENMPI/bin/mpif90
