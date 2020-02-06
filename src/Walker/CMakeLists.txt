cmake_minimum_required(VERSION 2.8.5)

project(Walker CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Walker
            Distributor.cpp
            Collector.cpp
            Integrator.cpp)

target_include_directories(Walker PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Main
                           ${QUINOA_SOURCE_DIR}/DiffEq
                           ${QUINOA_SOURCE_DIR}/IO
                           ${QUINOA_SOURCE_DIR}/Statistics
                           ${PROJECT_BINARY_DIR}/../Base
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PROJECT_BINARY_DIR}/../IO
                           ${PROJECT_BINARY_DIR}/../Walker
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${MKL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

addCharmModule( "distributor" "Walker" )
addCharmModule( "collector" "Walker" )
addCharmModule( "integrator" "Walker" )

# Add extra dependency of Walker on walkerCharmModule. This is required as
# Walker refers to the main Charm++ proxy defined in the Charm++ module
# walker (in Main/Walker.cpp).
add_dependencies( "Walker" "walkerCharmModule" )

# Add extra dependency of Integrator charm module on ParticleWriter charm
# module. This is required so that particlewriter.decl.h and
# particlewriter.def.h are generated before Distributor including those.
add_dependencies( "integratorCharmModule" "particlewriterCharmModule" )

set_target_properties(Walker PROPERTIES LIBRARY_OUTPUT_NAME quinoa_walker)

INSTALL(TARGETS Walker
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
