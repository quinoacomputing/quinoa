#!/bin/bash -eux
# vim: filetype=sh:
################################################################################
# 
# \file      script/run_tests.sh
# \author    J. Bakosi
# \date      Fri 29 Apr 2016 09:25:59 AM MDT
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
#  Note that it should also be fine putting the commands in this script into
#  .drone.yml, however, for some unknown reason, only the first command was
#  executed, so we call this script, running multiple test suites, used for
#  automated testing.
#
################################################################################

# Query number of CPUs
case "$OSTYPE" in
  darwin*)  CPUS=`sysctl -n hw.ncpu`;;
  linux*)   CPUS=`cat /proc/cpuinfo | grep MHz | wc -l`;;
esac

# Run unit test suite
./charmrun +p$CPUS --allow-run-as-root Main/unittest -v

# Run regression test suite (skip stringent tests that would run very long)
ctest -j$CPUS -LE stringent
