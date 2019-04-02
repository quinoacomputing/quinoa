cmake_minimum_required(VERSION 2.8.5)

project(PDE CXX)

add_library(PDE
            PDEStack.C
            Integrate/Quadrature.C
            Integrate/Initialize.C
            Integrate/Mass.C
            Integrate/Surface.C
            Integrate/Boundary.C
            Integrate/Volume.C
            Integrate/Source.C
            Integrate/Basis.C
            Integrate/Riemann/RiemannFactory.C
            Transport/Problem/CylAdvect.C
            Transport/Problem/GaussHump.C
            Transport/Problem/ShearDiff.C
            Transport/Problem/SlotCyl.C
            CompFlow/Physics/CGNavierStokes.C
            CompFlow/Problem/NLEnergyGrowth.C
            CompFlow/Problem/RayleighTaylor.C
            CompFlow/Problem/SodShocktube.C
            CompFlow/Problem/RotatedSodShocktube.C
            CompFlow/Problem/SedovBlastwave.C
            CompFlow/Problem/TaylorGreen.C
            CompFlow/Problem/UserDefined.C
            CompFlow/Problem/VorticalFlow.C
            Limiter.C
            ConfigureTransport.C
            ConfigureCompFlow.C
            ConfigureMultiMat.C)

target_include_directories(PDE PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/PDE
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

set_target_properties(PDE PROPERTIES LIBRARY_OUTPUT_NAME quinoa_pde)

INSTALL(TARGETS PDE
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
