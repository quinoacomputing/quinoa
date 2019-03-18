#!/bin/bash -e
################################################################################
#
# \file      script/update_copyright.sh
# \brief     Switch copyright year in all files
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
################################################################################

# Suggested argument directories:
#
#   quinoa: doc docker tests script src
#   quinoa-tpl: docker
#   quinoa-cmake: .

find doc -type f -not -name update_copyright.sh -exec sed -i 's/\\copyright.*/\\copyright 2012-2015 J. Bakosi,\n             2016-2018 Los Alamos National Security, LLC.,\n             2019 Triad National Security, LLC.\n             All rights reserved. See the LICENSE file for details./' {} +
