#!/bin/bash -eu
# vim: filetype=sh:
################################################################################
# 
# \file      script/run_tests.sh
# \author    J. Bakosi
# \date      Mon 29 Feb 2016 01:15:41 PM MST
# \copyright 2012-2016, Jozsef Bakosi.
# \brief     Run multiple test suites as part of automated testing
# \details   Run multiple test suites as part of automated testing.
#
#   Arguments to bash:
#   -e: shell will exit if any statement returns a non-true return value
#   -u: shell will exit if we try to use an uninitialised variable
#
#  This script will try to run the test suites in whatever build directory it is
#  called in.
#
################################################################################

# Query number of CPUs
CPUS=`cat /proc/cpuinfo | grep MHz | wc -l`

# Run unit test suite
./charmrun +p$CPUS Main/unittest -v

# Run regression test suite
ctest -j$CPUS -R GeneralizedDirichlet -E Big
