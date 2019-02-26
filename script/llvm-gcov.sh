#!/bin/bash
###############################################################################
#
# \file      script/llvm-gcov.sh
# \copyright 2012-2015, J. Bakosi, 2016-2019, Los Alamos National Security, LLC.
# \brief     Helper script for calling llvm-cov
#
################################################################################
exec llvm-cov-3.8 gcov "$@"
