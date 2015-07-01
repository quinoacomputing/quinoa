#!/usr/bin/env sh
#
# This script performs the refactoring from the old MOOCHO to MOOCHO
# as a Trilinos package proper.  This script should be run from the base directory of
# any source tree that you want to change.
#
# To run this script, you must have the variable $MOOCHO_BASE_DIR set
# to the base MOOCHO directory so that
#
#    $MOOCHO_BASE_DIR/moocho/src/ConstrainedOptPack/refactoring
#
# is where this file is.  Trilinos/commonTools/refactoring must also
# be added to your path.
#

token-replace-list-r $MOOCHO_BASE_DIR/moocho/src/ConstrainedOptPack/refactoring/new-header-includes.20051209.txt
