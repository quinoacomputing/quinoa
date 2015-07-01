#!/usr/bin/env sh
#
# This script performs the refactoring from the old MOOCHO to MOOCHO
# as a Trilinos package proper.  This script should be run from the base directory of
# any source tree that you want to change.
#
# To run this script, you must have the variable $MOOCHO_BASE_DIR set
# to the base MOOCHO directory so that
#
#    $MOOCHO_BASE_DIR/moocho/refactoring
#
# is where this file is.  Trilinos/commonTools/refactoring must also
# be added to your path.
#

$MOOCHO_BASE_DIR/moocho/src/MoochoUtilities/refactoring/from-old-to-trilinos.20051206.sh
$MOOCHO_BASE_DIR/moocho/src/IterationPack/refactoring/from-old-to-trilinos.20051207.sh
$MOOCHO_BASE_DIR/moocho/src/RTOpPack/refactoring/from-old-to-trilinos.20051207.sh
$MOOCHO_BASE_DIR/moocho/src/DenseLinAlgPack/refactoring/from-old-to-trilinos.20051208.sh
$MOOCHO_BASE_DIR/moocho/src/AbstractLinAlgPack/refactoring/from-old-to-trilinos.20051208.sh
$MOOCHO_BASE_DIR/moocho/src/NLPInterfacePack/refactoring/from-old-to-trilinos.20051209.sh
$MOOCHO_BASE_DIR/moocho/src/ConstrainedOptPack/refactoring/from-old-to-trilinos.20051209.sh
$MOOCHO_BASE_DIR/moocho/src/MoochoPack/refactoring/from-old-to-trilinos.20051209.sh
