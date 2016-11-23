#!/bin/bash -eux
# vim: filetype=sh:
################################################################################
# 
# \file      script/run_tests.sh
# \author    J. Bakosi
# \date      Wed 23 Nov 2016 11:14:18 AM MST
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Run multiple test suites as part of automated testing
# \details   Run multiple test suites as part of automated testing.
#
#   Arguments to bash:
#   -e: shell will exit if any statement returns a non-true return value
#   -u: shell will exit if we try to use an uninitialized variable
#   -x: shell will print each command to stdout before executing it
#
#  This script will try to run the test suites in whatever build directory it is
#  called in.
#
#  Arguments to script:
#  * No argument: Grab all available CPUs
#  * 1st optional argument: Use that as the number of CPUs
#  * 2nd optional argument: Forward argument to charmrun, e.g., mpirun
#
################################################################################

# If an argument is given, use that as the number of CPUs, if not grab them all
if [ "$#" -ge 1 ]; then
  CPUS=$1
else
  # Query number of CPUs
  case "$OSTYPE" in
    darwin*)  CPUS=`sysctl -n hw.ncpu`;;
    linux*)   CPUS=`cat /proc/cpuinfo | grep MHz | wc -l`;;
  esac
fi
echo "Will use $CPUS CPUs"

if [ "$#" -eq 2 ]; then
  CHARMRUN_ARG=$2
  echo "Will pass '$CHARMRUN_ARG' to charmrun"
else
  CHARMRUN_ARG=''
fi

# Configure parallel job runner
if [ ! -z ${NERSC_HOST:-} ]; then
  RUNNER=srun
  NCPUS_ARG=-n
else
  RUNNER=./charmrun
  NCPUS_ARG=+p
fi
# Run unit test suite
${RUNNER} ${NCPUS_ARG} $CPUS $CHARMRUN_ARG $PWD/Main/unittest -v

# Run regression test suite (skip stringent tests that would run very long)
ctest -j$CPUS -LE stringent
