cmake_minimum_required(VERSION 2.8.5)

project(Walker CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Walker
            Distributor.C
            Collector.C
            Integrator.C
)

addCharmModule( "distributor" "Walker" )
addCharmModule( "collector" "Walker" )
addCharmModule( "integrator" "Walker" )

# Add extra dependency of Walker on walkerCharmModule. This is required as
# Walker refers to the main Charm++ proxy defined in the Charm++ module
# walker (in Main/Walker.C).
add_dependencies( "Walker" "walkerCharmModule" )

set_target_properties(Walker PROPERTIES LIBRARY_OUTPUT_NAME quinoa_walker)

INSTALL(TARGETS Walker
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
