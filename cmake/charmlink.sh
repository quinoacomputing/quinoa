#!/bin/sh
comp=$1                                     # save charmc path and wrapper
shift 2                                     # skip the linker
for var in "$@"; do comp="$comp $var"; done # add rest of arguments
$comp                                       # link
