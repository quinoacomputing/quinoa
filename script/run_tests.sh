#!/bin/bash -eu
# vim: filetype=sh:
################################################################################
# 
# \file      script/run_tests.sh
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Run multiple test suites as part of automated testing
#
#   Arguments to bash:
#   -e: shell will exit if any statement returns a non-true return value
#   -u: shell will exit if we try to use an uninitialized variable
#   -x: shell will print each command to stdout before executing it
#
#   Run unit-, and regression test suites
#
#   run_tests.sh  [ncpus] [runner] [runner_ncpus_arg] [runner_args] [postfix_runner_args]
#
#   This script will try to run the test suites in whatever build directory it is
#   called in.
#
#   All arguments are optional. If no arguments are given, the script detects
#   the number of CPUs available and uses './charmrun' as the runner.
#
#   Arguments:
#
#     [ncpus] - Optionally specify the number of CPUs to use. Default: use all.
#
#     [runner] - Optionally specify the runner to prefix the unit test suite
#     with. Default: './charmrun'. Note that if a runner is specified the
#     argument [runner_ncpus_arg] must also be specified.
#
#     [runner_ncpus_arg] - Specify the argument used to specify the number of
#     CPUs for the runner, given by the [runner] argument.
#
#     [runner_args] - Extra runner arguments (before the exeuctable).
#
#     [postfix_runner_args] - Extra runner arguments (after the exeuctable).
#
#   Examples:
#
#     $ run_tests.sh
#     Detect number of CPUs and use ./charmrun to use them all.
#
#     $ run_tests.sh 12
#     Use 12 CPUs and use ./charmrun.
#
#     $ run_tests.sh 12 mpirun -n
#     Use 12 CPUs and use mpirun with -n specifying the number of CPUs.
#
#     $ run_tests.sh 12 mpirun -n "--bind-to none -oversubscribe"
#     Use 12 CPUs and use mpirun with -n specifying the number of CPUs and also
#     pass "--bind-to none -oversubscribe" to the runner.
#       - Note that the quotes are important so that multiple arguments are
#       interpreted as (prefix) runner arguments.
#
#     $ run_tests.sh 34 ./charmrun +p "--bind-to none -oversubscribe" "+ppn 17"
#     Use 36 CPUs distributed as 2 logical (compute) nodes x 17 worker threads
#     + 1 communication thread (2x18) in Charm++'s SMP mode using charmrun as
#     the runner with +p specifying the total number of worker threads and also
#     pass "--bind-to none -oversubscribe" to the runner.
#       - Note that the Charm++ SMP mode argument +ppn must be specified as a
#       postfix_runner_arg.
#       - Note that the number of worker threads must be a multiple of the
#       number of threads specified by the ppn argument.
#       - See also cmake/add_regression_test.cmake.
#
#  Note that only the [ncpus] argument affects the regression tests. The runner
#  and its optional extra arguments to the regression tests are configured by
#  cmake's RUNNER, RUNNER_NCPUS_ARG, and RUNNER_ARGS arguments. Example:
#
#  $ cmake -DRUNNER=mpirun -DRUNNER_NCPUS_ARG=-n -DRUNNER_ARGS="--bind-to none -oversubscribe"
#
################################################################################

# Error checking on the number of arguments
if [ "$#" -eq 2 ]; then
  echo -e "Error: 2 args are not allowed: if the 2nd [runner] arg is given, it must be followed by the 3rd [runner_ncpus_arg] arg.  See the source of this file for details."
  exit -1
fi

# If the 1st argument is given, use that as the number of CPUs, if not, detect
if [ "$#" -ge 1 ]; then
  NUMPES=$1
else
  case "$OSTYPE" in
    darwin*)  NUMPES=`sysctl -n hw.ncpu`;;
    linux*)   NUMPES=`cat /proc/cpuinfo | grep MHz | wc -l`;;
  esac
fi

# If the 2nd argument is given, there must be a 3rd one as well, setting the
# runner and its argument that sets the number of CPUs
if [ "$#" -ge 2 ]; then
  RUNNER=$2
  RUNNER_NCPUS_ARG=$3
else
  RUNNER=./charmrun
  RUNNER_NCPUS_ARG=+p
fi

# If the 4th argument is given, pass that to the runner
if [ "$#" -ge 4 ]; then
  RUNNER_ARGS=$4
else
  RUNNER_ARGS=''
fi

# If the 5th argument is given, pass that to the runner after the unittest executable
if [ "$#" -ge 4 ]; then
  POSTFIX_RUNNER_ARGS=$5
else
  POSTFIX_RUNNER_ARGS=''
fi

HARDWARE_NUMPES=${NUMPES}

echo "RUNNER: ${RUNNER}"
echo "RUNNER_NCPUS_ARG: ${RUNNER_NCPUS_ARG}"
echo "NUMPES: ${NUMPES}"
echo "RUNNER_ARGS: ${RUNNER_ARGS}"
echo "POSTFIX_RUNNER_ARGS: ${POSTFIX_RUNNER_ARGS}"

# Override HARDWARE_NUMPES with the number of CPUs used in hardware, used for
# ctest. In Charm++'s SMP mode NUMPES (used in ./charmrun's +p argument) is NOT
# the total number of PEs used in hardware by unittest because NUMPES does not
# contain the communication threads (by default one per compute node). Thus if
# +ppn is given in POSTFIX_RUNNER_ARGS, we parse out the ppn value and add the
# number of communication threads to NUMPES so that ctest uses the same number
# of PEs as unittest. For more on SMP, see also cmake/add_regression_test.cmake.
if [ -n "${POSTFIX_RUNNER_ARGS}" ]; then
  case "$POSTFIX_RUNNER_ARGS" in
    *ppn* ) PPN=`echo $POSTFIX_RUNNER_ARGS | sed -n -e 's/^.*ppn\s//p'`
            # If the runner is mpirun, -n specifies the number of compute nodes
            if [[ $RUNNER == *"mpirun"* ]]; then
              NUMNODES=${NUMPES}
            else
              NUMNODES=$(( ${NUMPES} / ${PPN} ))
            fi
            HARDWARE_NUMPES=$(( ${NUMNODES} * (${PPN}+1) ))
            echo "PPN: ${PPN}"
            echo "NUMNODES: ${NUMNODES}"
            echo "HARDWARE_NUMPES: ${HARDWARE_NUMPES}" ;;
  esac
fi

# Echo commands as they are executed starting from here
set -x

# Run unit test suite
${RUNNER} ${RUNNER_NCPUS_ARG} ${NUMPES} ${RUNNER_ARGS} $PWD/Main/unittest -v -q ${POSTFIX_RUNNER_ARGS}

# Run regression test suite (skip 'extreme' tests that would run very long)
ctest -j$HARDWARE_NUMPES --output-on-failure -LE extreme
