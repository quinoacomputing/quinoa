cmake_minimum_required(VERSION 2.8.5)

project(Walker CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Walker
            Distributor.C
            Integrator.C
)

addCharmModule( "distributor" "Walker" )
addCharmModule( "integrator" "Walker" )
removeWarnings( "Distributor.C;Integrator.C" )

# Add extra dependency of Walker on walkerCharmModule. This is required as
# Walker refers to the main Charm++ proxy defined in the Charm++ module
# walker (in Main/Walker.C).
add_dependencies( "Walker" "walkerCharmModule" )

INSTALL(TARGETS Walker
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
