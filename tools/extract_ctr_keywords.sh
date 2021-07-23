#!/bin/bash -eu
# -e: Exit immediately if a command exits with a non-zero status
# -u: Treat unset variables as an error when substituting
################################################################################
#
# \file      tools/extract_ctr_keywords.sh
# \brief     Generate documentation pages for control file keywords
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \details   This script runs an executable and extracts the help for each of
# the executable's control file keyword, describing their documentation,
# generating a documentation page with this information.
#
# A single command line argument is required: the name of the executable.
################################################################################

export exe=$1

if ( ../Main/$exe -h | grep '\-C' > /dev/null ); then

  # generate summary of all control file keywords
  export screen_out=$(../Main/$exe -C | awk "/$exe Control File Keywords/,0" | grep "^[[:space:]]\|^-" | sed 's/^/        /')
  perl -i -p0e 's/(\@section $ENV{exe}_ctr_list List of all control file keywords).*(\@section $ENV{exe}_ctr_detail)/$1\n\n$ENV{screen_out}\n\n$2/s' pages/${exe}_ctr.dox

  # generate detailed description on all keywords
  export keywords=$(../Main/$exe -C | awk "/$exe Control File Keywords/,0" | grep "^[[:space:]]" | awk '{print $1}')
  export detail=$(for k in $keywords; do l=$(echo $k | sed 's/+/plus/' | sed 's/-/minus/' | sed 's/\./_/'); ../Main/$exe -H $k ++quiet | grep -v '^Quinoa>' | sed "s/$exe control file keyword/@subsection ${exe}_ctr_kw_$l Keyword/" | sed 's/Expected type:/_Expected type:_/' | sed 's/Lower bound:/_Lower bound:_/' | sed 's/Upper bound:/_Upper bound:_/' | sed 's/Expected valid choices:/_Expected valid choices:_/' | sed 's/</\\</g' | sed 's/>/\\>/g' ; done)
  perl -i -p0e 's/(\@section $ENV{exe}_ctr_detail Detailed description of control file keywords).*(\*\/)/$1\n$ENV{detail}\n\n$2/s' pages/${exe}_ctr.dox

fi
