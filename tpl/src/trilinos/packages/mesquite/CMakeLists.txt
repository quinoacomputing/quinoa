IF (NOT DEFINED PROJECT_NAME)

  CMAKE_MINIMUM_REQUIRED( VERSION 2.4 )
  PROJECT(Mesquite)
  SET( PACKAGE_NAME "Mesquite" )

  MACRO(TRIBITS_ADD_OPTION_AND_DEFINE USER_OPTION_NAME MACRO_DEFINE_NAME DOCSTRING DEFAULT_VALUE)
    #MESSAGE("TRIBITS_ADD_OPTION_AND_DEFINE: '${USER_OPTION_NAME}' '${MACRO_DEFINE_NAME}' '${DEFAULT_VALUE}'")
    SET( ${USER_OPTION_NAME} "${DEFAULT_VALUE}" CACHE BOOL "${DOCSTRING}" )
    IF(NOT ${MACRO_DEFINE_NAME} STREQUAL "")
      SET(${MACRO_DEFINE_NAME} "${DEFAULT_VALUE}")
    ENDIF()
  ENDMACRO()

  SET( ${PACKAGE_NAME}_ENABLE_TESTS "OFF" CACHE BOOL "Build Tests" )
  IF (${PACKAGE_NAME}_ENABLE_TESTS)
    ENABLE_TESTING()
  ENDIF()

ELSE()

  INCLUDE(TribitsAddOptionAndDefine)
  INCLUDE(TribitsPackageMacros)

  #
  # A) Define the package
  #

  TRIBITS_PACKAGE(Mesquite DISABLE_STRONG_WARNINGS)
  SET( ${PACKAGE_NAME}_ENABLE_TESTS ${Trilinos_ENABLE_TESTS} )

ENDIF()


INCLUDE(CheckFunctionExists)
INCLUDE(CheckIncludeFile)
SET ( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake" )
INCLUDE(AutoconfHeader)

MACRO(APPEND_DEFINE DEFINE_LIST DEFINE_NAME)
  IF(${DEFINE_NAME})
    SET( ${DEFINE_LIST} "${${DEFINE_LIST}} -D${DEFINE_NAME}" )
  ENDIF()
ENDMACRO()

#
# B) Set up package-specific options
#

# Read the package version number specified as the second argument
# to the AC_INIT macro in a GNU Autoconf configure.in file.
#
# Define the following variables:
# VERSION_STRING:  The second argument to AC_INIT
# MAJOR_VERSION:   For a version string of the form m.n.p..., m
# MINOR_VERSION:   For a version string of the form m.n.p..., n
# PATCH_VERSION:   For a version string of the form m.n.p..., p

file( READ "configure.ac" configure_IN )
string( REGEX REPLACE "^.*AC_INIT *\\([^,]+, *([^, )]+).*$" "\\1" VERSION_STRING "${configure_IN}" )
string( REGEX REPLACE "^([0-9]+)(\\..*)?$" "\\1" MAJOR_VERSION "${VERSION_STRING}" )
string( REGEX REPLACE "^[0-9]+\\.([0-9]+)(\\..*)?$" "\\1" MINOR_VERSION "${VERSION_STRING}" )
if ( VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+" )
  string( REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" PATCH_VERSION "${VERSION_STRING}" )
endif ( VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+" )


SET(MSQ_VERSION ${VERSION_STRING})
SET(MSQ_VERSION_MAJOR ${MAJOR_VERSION})
SET(MSQ_VERSION_MINOR ${MINOR_VERSION})
SET(MSQ_VERSION_PATCH ${PATCH_VERSION})
SET(MSQ_VERSION_STRING "\"Mesquite ${MSQ_VERSION}\"")
SET(Mesquite_NAMESPACE "Mesquite${MAJOR_VERSION}" CACHE STRING "Mesquite namespace" )
SET(MESQUITE_NS ${Mesquite_NAMESPACE})
if (NOT MESQUITE_NS STREQUAL "Mesquite" )
  SET(MESQUITE_NS_ALIAS 1)
endif()

#check for specific functions
CHECK_FUNCTION_EXISTS(times     HAVE_TIMES)
CHECK_FUNCTION_EXISTS(clock     HAVE_CLOCK)
CHECK_FUNCTION_EXISTS(cbrt      HAVE_CBRT)
CHECK_FUNCTION_EXISTS(vsnprintf HAVE_VSNPRINTF)
CHECK_FUNCTION_EXISTS(vsprintf  HAVE_VSPRINTF)
CHECK_FUNCTION_EXISTS(finite    HAVE_FINITE)
CHECK_FUNCTION_EXISTS(_Isfinite HAVE__ISFINITE)

#check for header files
CHECK_INCLUDE_FILE(ieeefp.h    MSQ_HAVE_IEEEFP_H)
CHECK_INCLUDE_FILE(ieeefp.h    MSQ_HAVE_STDINT_H)
CHECK_INCLUDE_FILE(inttypes.h  MSQ_HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE(termios.h   HAVE_TERMIOS_H)
CHECK_INCLUDE_FILE(sys/ioctl.h HAVE_SYS_IOCTL_H)
CHECK_INCLUDE_FILE(dlfcn.h     HAVE_DLFCN_H)
CHECK_INCLUDE_FILE(fcntl.h     HAVE_FCNTL_H)

#check for ITAPS files
FIND_FILE( IMESH_DEFS  iMesh-Defs.inc ENV IMESH_DEFS  )
FIND_FILE( IGEOM_DEFS  iGeom-Defs.inc ENV IGEOM_DEFS  )
FIND_FILE( IREL_DEFS   iMesh-Defs.inc ENV IREL_DEFS   )
IF (TPL_ENABLE_MPI)
  FIND_FILE( IMESHP_DEFS iMesh-Defs.inc ENV IMESHP_DEFS )
ENDIF()

IF(NOT HAVE_TIMES AND NOT HAVE_CLOCK)
  MESSAGE( FATAL_ERROR
    "\n***"
    "\n*** ERROR: Neither the function times nor the function clock could be found. At least one of these needs to exist or Mesquite cannot build properly."
    "\n***\n"
    )
ENDIF()

IF (${PROJECT_NAME} STREQUAL "Trilinos")
  SET(Mesquite_ENABLE_SHARED ${TRILINOS_BUILD_SHARED_LIBS})
ELSE()
  TRIBITS_ADD_OPTION_AND_DEFINE(Mesquite_ENABLE_SHARED
    HAVE_Mesquite_SHARED
    "Build as shared library (e.g. Linux *.so, or Windows *.dll)."
    OFF )
ENDIF()

#this needs to be done more like the mesquite version where 
#there is a loop to determine the value. Not sure how to do
#this though. Likely through a try compile.
SET(${PACKAGE_NAME}_MSQ_FUNCTION "" CACHE STRING "Set c++ preprocessor macro for function name.")
SET(MSQ_FUNCTION ${${PACKAGE_NAME}_MSQ_FUNCTION})

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TRAP_FPE
  MSQ_TRAP_FPE
  "Enable to turn on floating point exception trapping."
  ON )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_DISABLE_DEBUG_ASSERTS
  NDEBUG
  "Disable asserts in debug mode."
  ON )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_USE_FUNCTION_TIMERS
  MSQ_USE_FUNCTION_TIMERS
  "Enable function timers."
  OFF )

IF(${PACKAGE_NAME}_ENABLE_DEBUG)
  SET(MSQ_ENABLE_DEBUG)
ENDIF()


#
# ) Output config data
#

SET( MESQUITE_DEFINES "" )
APPEND_DEFINE( MESQUITE_DEFINES MSQ_TRAP_FPE )
APPEND_DEFINE( MESQUITE_DEFINES MSQ_ENABLE_DEBUG )
APPEND_DEFINE( MESQUITE_DEFINES NDEBUG )
APPEND_DEFINE( MESQUITE_DEFINES MSQ_USE_FUNCTION_TIMERS )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_TIMES )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_CLOCK )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_CBRT )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_VSNPRINTF )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_VSPRINTF )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_FINITE )
APPEND_DEFINE( MESQUITE_DEFINES HAVE__ISFINITE )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_IEEEFP_H )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_INTTYPES_H )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_TERMIOS_H )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_SYS_IOCTL_H )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_DLFCN_H )
APPEND_DEFINE( MESQUITE_DEFINES HAVE_FCNTL_H )

IF (${PACKAGE_NAME}_ENABLE_ExodusII)
  SET(HAVE_LIBEXOIIV2C ON)
  SET( MSQ_USING_EXODUS MSQ_USING_EXODUS )
  APPEND_DEFINE( MESQUITE_DEFINES MSQ_USING_EXODUS )
ENDIF()



#
# C) Add the libraries, tests, and examples
#


# Add the list of include paths

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/testSuite)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/itaps)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/itaps)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/src/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/Control)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/MappingFunction)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/MappingFunction/Linear)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/MappingFunction/Lagrange)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/Mesh)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/Misc)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/ObjectiveFunction)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityAssessor)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityImprover)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityImprover/OptSolvers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityImprover/Relaxation)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Debug)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Shape)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Smoothness)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Volume)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Untangle)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/TMP)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/TargetCalculator)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/Misc)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/Shape)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/ShapeOrient)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/ShapeSize)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/ShapeSizeOrient)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/Size)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/Untangle)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/Wrappers)

ADD_SUBDIRECTORY(itaps)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(utils)

IF(NOT WIN32 AND NOT WIN64)
TRIBITS_EXCLUDE_FILES(
  meshFiles/2D/CUBIT/bad_circle.cbt
  meshFiles/2D/CUBIT/bad_circle.g
  meshFiles/2D/CUBIT/circle.g
  meshFiles/2D/CUBIT/cube-clip-corner.exo
  meshFiles/2D/CUBIT/hexagon.g
  meshFiles/2D/CUBIT/MesquiteGeo.g
  meshFiles/2D/vtk/quads/untangled/copy-morph-for-ricoh.cub.vtk
  meshFiles/2D/vtk/quads/untangled/four_more_quads.vtk
  meshFiles/2D/vtk/quads/untangled/halfCirc.vtk
  meshFiles/2D/vtk/quads/untangled/hole_scheme_opp_bias.vtk
  meshFiles/2D/vtk/quads/untangled/hole_scheme.vtk
  meshFiles/2D/vtk/quads/untangled/horse20x10_2_4.vtk
  meshFiles/2D/vtk/tris/tangled/horseshoe.vtk
  meshFiles/2D/vtk/mixed/untangled/hybrid_3quad_1tri_tangled.vtk
  meshFiles/2D/vtk/mixed/untangled/overt_hyb_1.vtk
  meshFiles/2D/vtk/mixed/untangled/overt_hyb_2.vtk
  meshFiles/2D/vtk/quads/untangled/polarMesh.vtk
  meshFiles/2D/vtk/quads/untangled/quad_2_not_xy.vtk
  meshFiles/2D/vtk/quads/tangled/rotsq.vtk
  meshFiles/2D/vtk/tris/untangled//shashkov_opt.vtk
  meshFiles/2D/vtk/quads/untangled/shashkov_quad_opt.vtk
  meshFiles/2D/vtk/tris/untangled//shashkov.vtk
  meshFiles/2D/vtk/tris/untangled/square_tri_2.vtk
  meshFiles/2D/vtk/quads/untangled/subdiv4x4_3_2_3x3.vtk
  meshFiles/2D/vtk/quads/untangled/tfi_chevron24x24.vtk
  meshFiles/2D/vtk/quads/untangled/tfi_horse10x4-14.vtk
  meshFiles/2D/vtk/tris/untangled/tri_5_not_xy.vtk
  meshFiles/2D/vtk/quads/untangled/unstrquad.vtk
  meshFiles/2D/vtk/vtk_to_gnu.cc
  meshFiles/2D/vtk/quads/untangled/ybiasqr24x24.vtk
  meshFiles/2D/vtk/N-Polygonal/poly2.vtk
  meshFiles/2D/vtk/N-Polygonal/poly3.vtk
  meshFiles/2D/vtk/quadratic/homogeneousPart.vtk
  meshFiles/2D/vtk/quadratic/horseShoe0.vtk
  meshFiles/2D/vtk/quadratic/horseShoe1.vtk
  meshFiles/2D/vtk/quadratic/quadPart.vtk
  meshFiles/2D/vtk/quadratic/triPart.vtk
  meshFiles/2D/vtk/quadratic/triPatch_0_0.vtk
  meshFiles/2D/vtk/quadratic/triPatch_0_7.vtk
  meshFiles/2D/vtk/quadratic/triPatch_0_8.vtk
  meshFiles/2D/vtk/quadratic/triPatch_1_0.vtk
  meshFiles/2D/vtk/quadratic/triPatch_1_1.vtk
  meshFiles/2D/vtk/quadratic/triPatch_1_2.vtk
  meshFiles/3D/CUBIT/brick_with_tet.cbt
  meshFiles/3D/CUBIT/cube_tets_2.tet
  meshFiles/3D/CUBIT/equilateral_tet_4.tet
  meshFiles/3D/vtk/pyramids/untangled/6-pyramid-brick.vtk
  meshFiles/3D/vtk/tets/untangled/cube-clip-12.vtk
  meshFiles/3D/vtk/tets/untangled/cube-clip-5.vtk
  meshFiles/3D/vtk/hexes/untangled/cube_hex_2.vtk
  meshFiles/3D/vtk/hexes/untangled/cube_hex_im.vtk
  meshFiles/3D/vtk/hexes/untangled/cube_hex_rm.vtk
  meshFiles/3D/vtk/hexes/untangled/cub_fru_scaled_100.vtk
  meshFiles/3D/vtk/hexes/untangled/cub_fru.vtk
  meshFiles/3D/vtk/hexes/untangled/cub_wedge_like.vtk
  meshFiles/3D/vtk/hexes/untangled/hex_2_with_bound.vtk
  meshFiles/3D/vtk/hexes/tangled/hex_5_tangled.vtk
  meshFiles/3D/vtk/hexes/tangled/pave1.e.vtk
  meshFiles/3D/vtk/hexes/tangled/pave2.e.vtk
  meshFiles/3D/vtk/hexes/tangled/pave3.e.vtk
  meshFiles/3D/vtk/hexes/tangled/pave_cyl.e.vtk
  meshFiles/3D/vtk/to-do/quad_sphere_ref.vtk
  meshFiles/3D/vtk/to-do/quad_sphere.vtk
  meshFiles/3D/vtk/hexes/tangled/roshan.vtk
  meshFiles/3D/vtk/hexes/untangled/single_hex.vtk
  meshFiles/3D/vtk/hexes/untangled/small_box_hex_1000.vtk
  meshFiles/3D/vtk/tets/untangled/tet_2_with_bound.vtk
  meshFiles/3D/vtk/tets/untangled/tet_4_tangled.vtk
  meshFiles/3D/vtk/tets/untangled/two_tets_shape.vtk
  meshFiles/3D/vtk/tets/tangled/two_tets.vtk
  meshFiles/3D/vtk/quadratic/cutCube.vtk
  meshFiles/3D/vtk/quadratic/hex3Dpart.vtk
  meshFiles/3D/vtk/quadratic/hexCyl.vtk
  meshFiles/3D/vtk/quadratic/sphereCube.vtk
  meshFiles/3D/vtk/quadratic/sphereCylinder_1194_inv.vtk
  meshFiles/3D/vtk/quadratic/sphereCylinder_22492_inv.vtk
  meshFiles/3D/vtk/quadratic/sphereCylinder_22492.vtk
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_TESTS)
  IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")

    TRIBITS_ADD_TEST_DIRECTORIES(testSuite)

    #TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

  ELSE()

    ADD_SUBDIRECTORY(testSuite)

  ENDIF()
ENDIF()

IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")
  #
  # D) Do standard postprocessing
  #

  TRIBITS_PACKAGE_POSTPROCESS()
ENDIF()
