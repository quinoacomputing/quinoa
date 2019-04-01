cmake_minimum_required(VERSION 2.8.5)

project(DiffEq CXX)

add_library(DiffEq
            DiffEqStack.C
)

target_include_directories(DiffEq PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/RNG
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Statistics
                           ${QUINOA_SOURCE_DIR}/DiffEq
                           ${QUINOA_SOURCE_DIR}/DiffEq/Beta
                           ${QUINOA_SOURCE_DIR}/DiffEq/Dirichlet
                           ${QUINOA_SOURCE_DIR}/DiffEq/WrightFisher
                           ${QUINOA_SOURCE_DIR}/DiffEq/OrnsteinUhlenbeck
                           ${QUINOA_SOURCE_DIR}/DiffEq/Gamma
                           ${QUINOA_SOURCE_DIR}/DiffEq/SkewNormal
                           ${QUINOA_SOURCE_DIR}/DiffEq/Velocity
                           ${QUINOA_SOURCE_DIR}/DiffEq/Position
                           ${QUINOA_SOURCE_DIR}/DiffEq/Dissipation
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${RNGSSE2_INCLUDE_DIRS}
                           ${MKL_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS})

set_target_properties(DiffEq PROPERTIES LIBRARY_OUTPUT_NAME quinoa_diffeq)

INSTALL(TARGETS DiffEq
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

# Configure build targets for stochastic differential equation types
add_subdirectory(Beta)
add_subdirectory(Dirichlet)
add_subdirectory(WrightFisher)
add_subdirectory(OrnsteinUhlenbeck)
add_subdirectory(Gamma)
add_subdirectory(SkewNormal)
add_subdirectory(Velocity)
add_subdirectory(Position)
add_subdirectory(Dissipation)
