#!/bin/bash
###############################################################################
#
# \file      script/llvm-gcov.sh
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Helper script for calling llvm-cov
#
################################################################################

exec llvm-cov-3.8 gcov "$@"
