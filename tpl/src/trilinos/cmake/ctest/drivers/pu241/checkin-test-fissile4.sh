#!/bin/bash

# Used to test Trilinos on any of the ORNL fissle 4 machines
#
# This script requires that modules be loaded by sourcing the script:
#
#    /opt/casl_vri_dev_env/fissile_four/build_scripts/load_official_dev_env.[sh,csh]
#
# You can source this script either in your shell startup script
# (e.g. .bash_profile) or you can source it manually whenever you need
# to set up to build VERA software.
#
# NOTE: This script does *NOT* automatically picks up any CASL VRI related
# extra repos and adds them to the checkin-test.py --extra-repos argument.  If
# you want to add extra repos, you can just pass in
# --extra-repos=Repo1,Repo2,...
#
# NOTE: This script automatically add SS extra builds so that one can
# simply run the script without having to specifiy extra builds and
# the right thing will happen.  Just make sure that you have cloned
# the right repos.
#
# NOTE: This script should not be directly modifed by typical CASL
# developers except, perhaps to add new extra builds.

EXTRA_ARGS=$@

# The default location for this directory tree is:
#
#  Trilinos.base
#    Trilinos    (your Trilinos soruce tree)
#    BUILDS
#      CHECKIN   (where you run this script from)
#
if [ "$TRILINOS_BASE_DIR" == "" ] ; then
  TRILINOS_BASE_DIR=../..
fi

TRILINOS_BASE_DIR_ABS=$(readlink -f $TRILINOS_BASE_DIR)

DRIVERS_BASE_DIR="$TRILINOS_BASE_DIR_ABS/Trilinos/cmake/ctest/drivers/pu241"

# Pakage in Trilinos to disable (mostly for auotmated CI server)
DISABLE_PACKAGES=PyTrilinos,TriKota,Optika,Sundance,Stokhos

# Check to make sure that the env has been loaded correctly
if [ "$CASL_VERA_OFFICIAL_DEV_ENV_LOADED" == "" ] ; then
  echo "Error, must source /opt/casl_vri_dev_env/fissile_four/build_scripts/load_official_dev_env.[sh,csh] before running checkin-test-fissile4.sh!"
  exit 1
fi

#
# Built-in Primary Stable (PS) builds (DO NOT MODIFY)
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.6.1-mpi-debug-ps-options.cmake'
" > MPI_DEBUG.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.6.1-serial-release-ps-options.cmake'
" > SERIAL_RELEASE.config

#
# Standard Secondary Stable (SS) builds (DO NOT MODIFY)
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.6.1-mpi-debug-ss-options.cmake'
" > MPI_DEBUG_SS.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.6.1-serial-release-ss-options.cmake'
" > SERIAL_RELEASE_SS.config

#
# Hybrid GCC 4.5.1 / Intel 11 extra builds
#
# Just pass these in with:
#
#    --extra-builds=MPI_DEBUG_451_11,SERIAL_RELEASE_451_11,MPI_DEBUG_SS_451_11,SERIAL_RELEASE_SS_451_11
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.5.1-mpi-debug-ps-options.cmake'
" > MPI_DEBUG_451_11.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.5.1-serial-release-ps-options.cmake'
" > SERIAL_RELEASE_451_11.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.5.1-mpi-debug-ss-options.cmake'
" > MPI_DEBUG_SS_451_11.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.5.1-serial-release-ss-options.cmake'
" > SERIAL_RELEASE_SS_451_11.config


#
# Extra intel builds added with --extra-builds=INTEL11_SERIAL_RELEASE,...
#
# NOTE: Since these are not required for pre-push testing purposes,
# CASL developers can modify these and push the changes.
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/intel-12.0.4-serial-release-ss-options.cmake'
-DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON
-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF
" > INTEL12_SERIAL_RELEASE.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/intel-12.0.4-serial-debug-ss-options.cmake'
-DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON
-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF
" > INTEL12_SERIAL_DEBUG.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/intel-11.1.064-serial-release-ss-options.cmake'
-DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON
-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF
" > INTEL11_SERIAL_RELEASE.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/intel-11.1.064-serial-debug-ss-options.cmake'
-DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON
-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF
" > INTEL11_SERIAL_DEBUG.config

#
# Invocation
#

$TRILINOS_BASE_DIR/Trilinos/checkin-test.py \
-j16 \
--ctest-timeout=180 \
--ss-extra-builds=MPI_DEBUG_SS,SERIAL_RELEASE_SS \
--disable-packages=$DISABLE_PACKAGES \
--skip-case-no-email \
$EXTRA_ARGS  


# NOTE: By default we use 16 processes which is 1/2 of the 32 processes on a
# fissile 4 machine.  This way two people can build and test without taxing
# the machine too much.

# Here is the template for excluding tests that can be aded to
# checkin-test.py call above
#
#   --ctest-options="-E '(Test1|Test2)'" \
