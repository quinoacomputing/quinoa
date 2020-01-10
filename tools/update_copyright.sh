#!/bin/bash -e
################################################################################
#
# \file      tools/update_copyright.sh
# \brief     Switch copyright year in all files
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
################################################################################

# Suggested argument directories:
#
#   quinoa: doc tests tools src cmake
#   quinoa-tpl: docker
#   quinoa-cmake: .

find doc -type f -not -name update_copyright.sh -exec sed -i 's/\\copyright.*/\\copyright 2012-2015 J. Bakosi,\n             2016-2018 Los Alamos National Security, LLC.,\n             2019 Triad National Security, LLC.\n             All rights reserved. See the LICENSE file for details./' {} +
