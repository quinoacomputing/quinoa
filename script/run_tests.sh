#!/bin/bash -eux
# vim: filetype=sh:
################################################################################
# 
# \file      script/run_tests.sh
# \author    J. Bakosi
# \date      Tue 08 Mar 2016 08:04:35 AM MST
# \copyright 2012-2016, Jozsef Bakosi.
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
#  Note that it should also be fine putting the commands in this script into
#  .drone.yml, however, for some unknown reason, only the first command was
#  executed, so we call this script, running multiple test suites, used ofr
#  automated testing.
#
################################################################################

# Query number of CPUs
CPUS=`cat /proc/cpuinfo | grep MHz | wc -l`

# Run unit test suite
./charmrun +p$CPUS Main/unittest -v

# Run regression test suite (skip stringent tests that would run very long)
ctest -j$CPUS -LE stringent
