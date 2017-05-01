#!/bin/bash -eux
# vim: filetype=sh:
################################################################################
# 
# \file      script/run_tests.sh
# \author    J. Bakosi
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
#  * 2nd optional argument: Forward argument to runner, e.g., mpirun
#
################################################################################

# If an argument is given, use that as the number of CPUs, if not grab them all
if [ "$#" -ge 1 ]; then
  NUMPES=$1
else
  # Query number of CPUs
  case "$OSTYPE" in
    darwin*)  NUMPES=`sysctl -n hw.ncpu`;;
    linux*)   NUMPES=`cat /proc/cpuinfo | grep MHz | wc -l`;;
  esac
fi

if [ "$#" -eq 2 ]; then
  RUNNER_ARGS=$2
  echo "Will pass '$RUNNER_ARGS' to runner"
else
  RUNNER_ARGS=''
fi

# Configure parallel job runner
if [ ! -z ${NERSC_HOST:-} ]; then
  RUNNER=srun
  RUNNER_NCPUS_ARG=-n
else
  RUNNER=./charmrun
  RUNNER_NCPUS_ARG=+p
fi
# Run unit test suite
${RUNNER} ${RUNNER_NCPUS_ARG} ${NUMPES} ${RUNNER_ARGS} $PWD/Main/unittest -v

# Run regression test suite (skip 'extreme' tests that would run very long)
ctest -j$NUMPES --output-on-failure -LE extreme
