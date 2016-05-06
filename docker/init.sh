#!/bin/bash
################################################################################
#
# \file      docker/init.sh
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \date      Fri 06 May 2016 06:46:10 AM MDT
# \brief     Setup environment modules (for use in automated testing)
#
################################################################################

set -e

source /usr/share/modules/init/bash

exec "$@"
