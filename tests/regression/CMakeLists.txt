################################################################################
#
# \file      tests/regression/CMakeLists.txt
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Cmake code common to all regression testing
#
################################################################################

# Set paths required for running tests
set(TEST_RUNNER ${CMAKE_SOURCE_DIR}/../cmake/test_runner.cmake)

# Find gmsh executable
find_package(Gmsh)
if(ENABLE_MESHCONV AND NOT GMSH_FOUND)
  message(WARNING "Gmsh not found, meshconv regression tests will not be rigorous")
endif()

# Find valgrind executable
set(ENABLE_VALGRIND false CACHE BOOL "Run all regression tests using valgrind")
if (ENABLE_VALGRIND)
  find_program(VALGRIND valgrind)
  if (VALGRIND)
    set(USE_VALGRIND true)
    message(STATUS "Will run all regression tests using valgrind")
  else()
    message(WARNING "Valgrind was enabled but NOT found")
  endif()
endif()

# Include function used to add regression tests
include(add_regression_test)

# Include regression test directories

# Walker regression tests
if (ENABLE_WALKER)
  message(STATUS "Adding regression tests for ${WALKER_EXECUTABLE}")
  add_subdirectory(walker/Dirichlet)
  add_subdirectory(walker/MixDirichlet)
  add_subdirectory(walker/GeneralizedDirichlet)
  add_subdirectory(walker/SkewNormal)
  add_subdirectory(walker/DiagOrnsteinUhlenbeck)
  add_subdirectory(walker/OrnsteinUhlenbeck)
  add_subdirectory(walker/Beta)
  add_subdirectory(walker/Gamma)
  add_subdirectory(walker/NumFracBeta)
  add_subdirectory(walker/MassFracBeta)
  add_subdirectory(walker/MixMassFracBeta)
  add_subdirectory(walker/Velocity)
endif()

# RNGTest regression tests
if(ENABLE_RNGTEST)
  message(STATUS "Adding regression tests for ${RNGTEST_EXECUTABLE}")
  add_subdirectory(rngtest)
endif()

# MeshConv regression tests
if(ENABLE_MESHCONV)
  message(STATUS "Adding regression tests for ${MESHCONV_EXECUTABLE}")
  add_subdirectory(meshconv/gmsh_output)
  add_subdirectory(meshconv/netgen_output)
  add_subdirectory(meshconv/exo_output)
endif()

# Inciter regression tests
if(ENABLE_INCITER)
  message(STATUS "Adding regression tests for ${INCITER_EXECUTABLE}")
  add_subdirectory(inciter/transport/SlotCyl)
  add_subdirectory(inciter/transport/SlotCyl/fct)
  add_subdirectory(inciter/transport/SlotCyl/cfl)
  add_subdirectory(inciter/transport/SlotCyl/asynclogic)
  add_subdirectory(inciter/transport/ShearDiff)
  add_subdirectory(inciter/transport/GaussHump)
  add_subdirectory(inciter/transport/CylAdvect)
  add_subdirectory(inciter/compflow/Euler/VorticalFlow)
  add_subdirectory(inciter/compflow/Euler/TaylorGreen)
  add_subdirectory(inciter/compflow/Euler/RayleighTaylor)
  add_subdirectory(inciter/compflow/Euler/NLEnergyGrowth)
  add_subdirectory(inciter/compflow/Euler/GaussHump)
  add_subdirectory(inciter/mesh_refinement/t0ref)
  add_subdirectory(inciter/mesh_refinement/dtref)
  add_subdirectory(inciter/compflow/Euler/SodShocktube)
  add_subdirectory(inciter/compflow/Euler/SedovBlastwave)
  add_subdirectory(inciter/multimat/InterfaceAdvection)
  add_subdirectory(inciter/multimat/SodShocktube)
  add_subdirectory(inciter/restart)
endif()
