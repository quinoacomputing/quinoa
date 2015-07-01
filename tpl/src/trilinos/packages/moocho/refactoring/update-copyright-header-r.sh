#!/bin/sh
find . -name "*.hpp" -exec $MOOCHO_BASE_DIR/moocho/refactoring/update-copyright-header.pl $MOOCHO_BASE_DIR/moocho/refactoring/c++-copyright-header.txt {} {} \;
find . -name "*.cpp" -exec $MOOCHO_BASE_DIR/moocho/refactoring/update-copyright-header.pl $MOOCHO_BASE_DIR/moocho/refactoring/c++-copyright-header.txt {} {} \;
find . -name "*.h" -exec $MOOCHO_BASE_DIR/moocho/refactoring/update-copyright-header.pl $MOOCHO_BASE_DIR/moocho/refactoring/c++-copyright-header.txt {} {} \;
find . -name "*.c" -exec $MOOCHO_BASE_DIR/moocho/refactoring/update-copyright-header.pl $MOOCHO_BASE_DIR/moocho/refactoring/c++-copyright-header.txt {} {} \;
find . -name "*.ud" -exec $MOOCHO_BASE_DIR/moocho/refactoring/update-copyright-header.pl $MOOCHO_BASE_DIR/moocho/refactoring/c++-copyright-header.txt {} {} \;
find . -name "Makefile*" -exec $MOOCHO_BASE_DIR/moocho/refactoring/update-copyright-header.pl $MOOCHO_BASE_DIR/moocho/refactoring/makefile-copyright-header.txt {} {} \;
find . -name "configure*" -exec $MOOCHO_BASE_DIR/moocho/refactoring/update-copyright-header.pl $MOOCHO_BASE_DIR/moocho/refactoring/makefile-copyright-header.txt {} {} \;
