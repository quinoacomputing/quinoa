#!/bin/bash

# Used to test Trilinos on any of the ORNL CASL Fissile/Spy machines
#
# This script requires that the VERA dev env be loaded by sourcing the script:
#
#  . /projects/vera/gcc-4.8.3/load_dev_env.[sh,csh]
#
# You can source this script either in your shell startup script
# (e.g. .bash_profile) or you can source it manually whenever you need to set
# up to build VERA software.
#
# NOTE: This script should not be directly modified by typical CASL
# developers except, perhaps to add new extra builds.
#
# NOTE: You can pass through arguments with spaces using quotes like:
#
#   --ctest-options="-E '(Test1|Test2)'"
#
# and it will preserve the spaces correctly.  If you want to pass along
# quotes, you have to escape them like:
#
#   --ctest-options="-E \"(Test1|Test2)\""
#
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

DRIVERS_BASE_DIR="$TRILINOS_BASE_DIR_ABS/Trilinos/cmake/ctest/drivers/fissile4"

# Packages in Trilinos to disable (mostly for auotmated CI server)
DISABLE_PACKAGES=PyTrilinos,Pliris,Claps,STK,TriKota

# Check to make sure that the env has been loaded correctly
if [ "$LOADED_TRIBITS_DEV_ENV" != "gcc-4.8.3" ] ; then
  echo "Error, must source /projects/vera/gcc-4.8.3/load_dev_env.[sh,csh] before running checkin-test-vera.sh!"
  exit 1
fi

echo "
-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
" > COMMON.config

#
# Built-in Primary Tested (PT) --default-builds (DO NOT MODIFY)
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.8.3-base-options.cmake,$DRIVERS_BASE_DIR/trilinos-tpls-gcc.4.8.3.cmake'
" > MPI_DEBUG.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.8.3-base-options.cmake,$DRIVERS_BASE_DIR/trilinos-tpls-gcc.4.8.3.cmake'
" > SERIAL_RELEASE.config

#
# Standard Secondary Tested (ST) --st-extra-builds (DO NOT MODIFY)
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.8.3-base-options.cmake,$DRIVERS_BASE_DIR/trilinos-tpls-gcc.4.8.3.cmake'
-DCMAKE_BUILD_TYPE=RELEASE
-DTrilinos_ENABLE_DEBUG=ON
-DTPL_ENABLE_MPI=ON
" > MPI_DEBUG_ST.config

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.8.3-base-options.cmake,$DRIVERS_BASE_DIR/trilinos-tpls-gcc.4.8.3.cmake'
-DCMAKE_BUILD_TYPE=RELEASE
-DTrilinos_ENABLE_DEBUG=OFF
-DTPL_ENABLE_MPI=OFF
" > SERIAL_RELEASE_ST.config

#
# --extra-builds
#

echo "
-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH='$DRIVERS_BASE_DIR/gcc-4.8.3-base-options.cmake,$DRIVERS_BASE_DIR/trilinos-tpls-gcc.4.8.3.cmake'
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTrilinos_ENABLE_DEBUG:BOOL=OFF
" > MPI_RELEASE.config

# Create local defaults file if one does not exist
_LOCAL_CHECKIN_TEST_DEFAULTS=local-checkin-test-defaults.py
if [ -f $_LOCAL_CHECKIN_TEST_DEFAULTS ] ; then
  echo "File $_LOCAL_CHECKIN_TEST_DEFAULTS already exists, leaving it!"
else
  echo "Creating default file $_LOCAL_CHECKIN_TEST_DEFAULTS!"
  echo "
defaults = [
  \"-j16\",
  \"--ctest-timeout=180\",
  \"--st-extra-builds=MPI_DEBUG_ST,SERIAL_RELEASE_ST\",
  \"--disable-packages=$DISABLE_PACKAGES\",
  \"--skip-case-no-email\",
  \"--ctest-options=-E '(MueLu_ParameterListInterpreterEpetra|MueLu_ParameterListInterpreterTpetra|Belos_pseudo_ptfqmr_hb_1_MPI_4|Belos_pseudo_ptfqmr_hb_3_MPI_4)'\",
  ]
  " > $_LOCAL_CHECKIN_TEST_DEFAULTS
fi

#
# Invocation
#

$TRILINOS_BASE_DIR/Trilinos/checkin-test.py \
"$@"


# --ctest-options="-E '(Piro_AnalysisDriver|Stokhos_Linear2D_Diffusion_GMRES_KLR|Panzer_STK_ResponseLibraryTest|MueLu_|Amesos2_|Rythmos_ImplicitRK_UnitTest_MPI_1|SEACASExodus_exodus_unit_tests|Intrepid_test_Discretization_Basis_HGRAD_TRI_Cn_FEM_Test_02_MPI_1|Intrepid_test_Discretization_Basis_HDIV_TET_In_FEM_Test_02_MPI_1|Intrepid_test_Discretization_Basis_HGRAD_TET_Cn_FEM_Test_02_MPI_1|Sundance_BesselTest2D_MPI_1|ThyraTpetraAdapters_TpetraThyraWrappersUnitTests_serial|Ifpack2_RILUKSingleProcessUnitTests)'"

# NOTE: By default we use 16 processes which is 1/2 of the 32 processes on a
# fissile 4 machine.  This way two people can build and test without taxing
# the machine too much.
