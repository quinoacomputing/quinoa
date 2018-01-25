################################################################################
#
# \file      script/update_copyright.sh
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Switch copyright year in all files
# \details   Suggested arguments:
#   quinoa: doc docker regression script src
#   quinoa-tpl: docker
#   quinoa-cmake: .
#
################################################################################

# Require at least one argument
die () { echo >&2 "$@"; exit 1; }
[ -z $1 ] && die "Usage: $0 <directories>"

find "$@" -type f -not -name update_copyright.sh -exec sed -i 's/opyright 2012-2015, J. Bakosi, 2016-2017/opyright 2012-2015, J. Bakosi, 2016-2018/g' {} +
find "$@" -type f -not -name update_copyright.sh -exec sed -i 's/Copyright (c) 2016-2017/Copyright (c) 2016-2018/g' {} +
find "$@" -type f -not -name update_copyright.sh -exec sed -i 's/Copyright 2016-2017/Copyright 2016-2018/g' {} +

find LICENSE -type f -exec sed -i 's/Copyright (c) 2016-2017/Copyright (c) 2016-2018/g' {} +
find LICENSE -type f -exec sed -i 's/Copyright 2016-2017/Copyright 2016-2018/g' {} +
