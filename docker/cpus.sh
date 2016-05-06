#!/bin/bash
################################################################################
#
# \file      docker/cpus.sh
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \date      Fri 06 May 2016 06:46:02 AM MDT
# \brief     Return the number of CPUs (for use in automated testing)
#
################################################################################

echo `cat /proc/cpuinfo | grep MHz | wc -l`
