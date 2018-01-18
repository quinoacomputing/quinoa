#!/bin/sh -ex
# ################################################################################
#
# \file      script/cobalt-runjob.sh
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Helper script for the cobalt scheduler
#
################################################################################
RUNNER_NCPUS_ARG=$1
NUMPES=$2
TEST_EXECUTABLE=$3
shift 3
TEST_EXECUTABLE_ARGS="$@"
runjob --block $COBALT_PARTNAME $RUNNER_NCPUS_ARG $NUMPES --ranks-per-node 16 --verbose 2 : $TEST_EXECUTABLE $TEST_EXECUTABLE_ARGS
