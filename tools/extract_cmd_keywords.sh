#!/bin/bash -eu
# -e: Exit immediately if a command exits with a non-zero status
# -u: Treat unset variables as an error when substituting
################################################################################
#
# \file      tools/extract_cmd_keywords.sh
# \brief     Generate documentation pages for command line keywords
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \details   This script runs an executable and extracts the help for each of
# the executable's command line keyword, describing their documentation,
# generating a documentation page with this information.
#
# A single command line argument is required: the name of the executable.
################################################################################

export exe=$1

# generate summary of all command line keywords
export screen_out=$(../Main/$exe -h | awk "/$exe Command-line Parameters/,0" | grep "^[[:space:]]\|^-" | sed 's/^/        /')
perl -i -p0e 's/(\@section $ENV{exe}_cmd_list List of all command line parameters).*(\@section $ENV{exe}_cmd_detail)/$1\n\n$ENV{screen_out}\n\n$2/s' pages/${exe}_cmd.dox

# generate detailed description on all keywords
export commands=$(../Main/$exe -h | awk "/$exe Command-line Parameters/,0" | grep "^[[:space:]]\|^-" | awk '{print $2}' | sed 's/\-\-//')
export detail=$(for c in $commands; do ../Main/$exe -H $c ++quiet | grep -v '^Quinoa>' | sed "s/$exe command-line keyword/@subsection ${exe}_cmd_kw_$c Keyword/" | sed 's/--/\\--/' | sed 's/Expected type:/_Expected type:_/' | sed 's/Lower bound:/_Lower bound:_/' | sed 's/Upper bound:/_Upper bound:_/' | sed 's/Expected valid choices:/_Expected valid choices:_/'; done)
perl -i -p0e 's/(\@section $ENV{exe}_cmd_detail Detailed description of command line parameters).*(\*\/)/$1\n$ENV{detail}\n\n$2/s' pages/${exe}_cmd.dox
