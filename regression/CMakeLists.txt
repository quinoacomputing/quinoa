################################################################################
#
# \file      src/regression/CMakeLists.txt
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Cmake code common to all regression testing
#
################################################################################

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Set paths required for running tests
set(TEST_RUNNER ${CMAKE_SOURCE_DIR}/../cmake/test_runner.cmake)

# Find gmsh executable
find_program(GMSH gmsh)
if(NOT GMSH)
  message(STATUS "Gmsh not found, meshconv regression tests will not be rigorous")
endif()

# Include function used to add regression tests
include(add_regression_test)

# Include regression test directories

# Walker regression tests
add_subdirectory(walker/Dirichlet)
add_subdirectory(walker/GeneralizedDirichlet)
add_subdirectory(walker/SkewNormal)
add_subdirectory(walker/DiagOrnsteinUhlenbeck)
add_subdirectory(walker/OrnsteinUhlenbeck)
add_subdirectory(walker/Beta)
add_subdirectory(walker/Gamma)
add_subdirectory(walker/NumFracBeta)
add_subdirectory(walker/MassFracBeta)
add_subdirectory(walker/MixMassFracBeta)

# RNGTest regression tests
if(HAS_TESTU01)
  add_subdirectory(rngtest)
endif()

# MeshConv regression tests
add_subdirectory(meshconv/gmsh_output)
add_subdirectory(meshconv/netgen_output)
add_subdirectory(meshconv/exo_output)

# Inciter regression tests
add_subdirectory(inciter/transport/SlotCyl/asynclogic)
add_subdirectory(inciter/transport/SlotCyl/fct)
add_subdirectory(inciter/transport/SlotCyl/cfl)
add_subdirectory(inciter/transport/ShearDiff)
add_subdirectory(inciter/transport/GaussHump)
add_subdirectory(inciter/compflow/Euler/VorticalFlow)
add_subdirectory(inciter/compflow/Euler/TaylorGreen)
add_subdirectory(inciter/compflow/Euler/RayleighTaylor)
add_subdirectory(inciter/compflow/Euler/NLEnergyGrowth)
#add_subdirectory(inciter/amr/initial/uniform)
