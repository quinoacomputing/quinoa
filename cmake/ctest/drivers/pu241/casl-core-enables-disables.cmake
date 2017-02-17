#
# This file is meant to be included in the configuration of Trilinos to
# disable packages and other critical varaibles internally.  These excludes
# can be overwritten on the command line in the CMake cache.
#
# These options are set for *all* compilers for all builds so they have to be
# 100% general!
#

# Must include first so that it defines ${PROJECT_NAME}_EXCLUDE_PACKAGES
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-exclude-trilinos-packages.cmake)

# Put in hard disables for excluded packages
FOREACH(EXCLUDED_PACKAGE ${${PROJECT_NAME}_EXCLUDE_PACKAGES})
  SET(${PROJECT_NAME}_ENABLE_${EXCLUDED_PACKAGE} OFF CACHE BOOL
    "Disabled in casl-core-enables-disables.cmake")
ENDFOREACH()

# Turn off float and complex testing because CASL does not need them
SET(Teuchos_ENABLE_FLOAT OFF CACHE BOOL "")
SET(Teuchos_ENABLE_COMPLEX OFF CACHE BOOL "")

# We don't want or need HDF5 support in EpetraExt
SET(EpetraExt_ENABLE_HDF5 OFF CACHE BOOL "Disabled in casl-core-enables-disables.cmake")

# Turn off STK tests since they are constantly failing.  NOTE: Since CASL is
# not developing on STK, only using it, this should not represent a big risk
# for STK or CASL.
SET(STK_ENABLE_TESTS OFF CACHE BOOL "Disabled in casl-core-enables-disables.cmake")
SET(STK_ENABLE_EXAMPLES OFF CACHE BOOL "Disabled in casl-core-enables-disables.cmake")

# Turn off this failing Rythmos test (see Trilinos bug 5485)
SET(Rythmos_ImplicitRK_ConvergenceTest_MPI_1_DISABLE TRUE)

# Turn off failing ML tests (see Trilinos bug 5537)
SET(ML_Blackboard_MPI_4_DISABLE TRUE)
SET(ML_AdaptiveSA_MPI_4_DISABLE TRUE)

# Turn off some failing Belos and Anasazi tests (see Trilinos bugs 5382 and 5383)
SET(Belos_Tpetra_MVOPTester_complex_test_MPI_4_DISABLE TRUE)
SET(Belos_Tpetra_MVOPTester_complex_test_DISABLE TRUE)
SET(Anasazi_Tpetra_MVOPTester_MPI_4_DISABLE TRUE)
SET(Anasazi_Tpetra_MVOPTester_DISABLE TRUE)

# Don't allow Optika GUI tests since they just seem to work and are not tested
# (see discussion in commit message).
SET(Optika_DO_GUI_UNIT_TESTS OFF CACHE BOOL "")

# Turn on configure timing
SET(${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING ON CACHE BOOL "")

# Don't create *Config.cmake files since they are massively expensive to create
SET(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES OFF CACHE BOOL "")
