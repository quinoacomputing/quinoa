

# dummy change to test secondary git repo
# another dummy change

#
# A) Define the package
#

TRIBITS_PACKAGE(Sundance  DISABLE_CIRCULAR_REF_DETECTION_FAILURE)

#
# B) Set up package-specific options
#

#IF (NOT (${Trilinos_DATA_DIR} STREQUAL "NOTFOUND"))
#   SET( HAVE_TRILINOS_DATA ON)
#ENDIF()
SET (HAVE_TRILINOS_DATA OFF)


TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Python
  HAVE_SUNDANCE_PYTHON
  "Enable python functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MOOCHO
  HAVE_SUNDANCE_MOOCHO
  "Enable Sundance/Moocho."
  ON )


IF(${PACKAGE_NAME}_ENABLE_Stokhos)
  SET(HAVE_SUNDANCE_STOKHOS_DEFAULT ON)
ELSE()
  SET(HAVE_SUNDANCE_STOKHOS_DEFAULT OFF)
ENDIF()
#Adding a define for when Stokhos is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_SUNDANCE_STOKHOS ${HAVE_SUNDANCE_STOKHOS_DEFAULT})

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_BROKEN_CODE
  HAVE_SUNDANCE_BROKEN_CODE
  "Enable broken code. For expert developer use only."
  OFF )

#Enable ExodusII reader in Sundance (via SEACAS package).
#Adding a define for when SEACASExodus is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_SUNDANCE_EXODUS "${${PACKAGE_NAME}_ENABLE_SEACASExodus}")

TRIBITS_ADD_OPTION_AND_DEFINE(Sundance_ENABLE_EXTERNAL_CHACO
  HAVE_SUNDANCE_CHACO
  "Enable file-based chaco partitioning in Sundance (requires chaco executable in path)."
  OFF )

# ========== Peano configuration ===============
IF(TPL_ENABLE_Peano)
  #Adding a define for when the TPL Peano is enabled.
  #This replaces an option that was reusing an existing name.
  SET(HAVE_SUNDANCE_PEANO "${${PACKAGE_NAME}_ENABLE_Peano}")

  INCLUDE_DIRECTORIES(${Peano_INCLUDE_DIRS}) 
  
  IF (TPL_ENABLE_Peano_NO_2D)
    TRIBITS_ADD_OPTION_AND_DEFINE(TPL_ENABLE_Peano_NO_2D
      HAVE_SUNDANCE_PEANO_NO_2D
     "Turn OFF Peano 2D"
     ON )
  ELSE()
    # LINK_DIRECTORIES( ${Peano_LIBRARY_DIRS}/libpeano2D.a )
    # PEANO_ADD_LIBRARY(peano2D)
  ENDIF()
  
  IF (TPL_ENABLE_Peano_NO_3D)
    TRIBITS_ADD_OPTION_AND_DEFINE(TPL_ENABLE_Peano_NO_3D
      HAVE_SUNDANCE_PEANO_NO_3D
     "Turn OFF Peano 3D"
     ON )
  ELSE()
    # LINK_DIRECTORIES( ${Peano_LIBRARY_DIRS}/libpeano3D.a )
    # PEANO_ADD_LIBRARY(peano3D) 
  ENDIF()
     
ENDIF()
# ============ end of Peano configuration ================

IF(${Proj}_VERBOSE_CONFIGURE)
  PRINT_VAR(HAVE_SUNDANCE_EXODUS)
  PRINT_VAR(HAVE_SUNDANCE_PYTHON)
  PRINT_VAR(HAVE_SUNDANCE_CHACO)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

#MESSAGE("adding Playa")
ADD_SUBDIRECTORY(Playa)

#MESSAGE("adding src-utils")
ADD_SUBDIRECTORY(src-utils)
#MESSAGE("adding test-utils")
TRIBITS_ADD_TEST_DIRECTORIES(tests-utils)

#ADD_SUBDIRECTORY(src-solvers)
#TRIBITS_ADD_TEST_DIRECTORIES(tests-solvers)

#MESSAGE("adding src-interop")
ADD_SUBDIRECTORY(src-interop)

#MESSAGE("adding src-core")
ADD_SUBDIRECTORY(src-core)
#MESSAGE("adding tests-core")
TRIBITS_ADD_TEST_DIRECTORIES(tests-core)

#MESSAGE("adding src-std-mesh")
ADD_SUBDIRECTORY(src-std-mesh)
#MESSAGE("adding tests-std-mesh")
TRIBITS_ADD_TEST_DIRECTORIES(tests-std-mesh)

#MESSAGE("adding src-std-fwk")
ADD_SUBDIRECTORY(src-std-framework)
#MESSAGE("adding tests-std-fwk")
TRIBITS_ADD_TEST_DIRECTORIES(tests-std-framework)

#MESSAGE("adding src-pdeopt")
ADD_SUBDIRECTORY(src-pdeopt)
#MESSAGE("adding tests-pdeopt")
TRIBITS_ADD_TEST_DIRECTORIES(tests-pdeopt)

IF(Sundance_ENABLE_Python)
  MESSAGE("adding python")
  ADD_SUBDIRECTORY(python)
ENDIF()

#MESSAGE("adding ex-tut")
TRIBITS_ADD_TEST_DIRECTORIES(examples)


#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
