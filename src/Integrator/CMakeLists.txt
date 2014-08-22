cmake_minimum_required(VERSION 2.8.5)

project(Integrator CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Integrator
            Distributor.C
            Integrator.C
)

addCharmModule( "distributor" "Integrator" )
addCharmModule( "integrator" "Integrator" )
removeWarnings( "Distributor.C;Integrator.C")

# Add extra dependency of Integrator on quinoaCharmModule. This is required as
# Integrator refers to the main Charm++ proxy defined in the Charm++ module
# quinoa (in Main/Quinoa.C).
add_dependencies( "Integrator" "quinoaCharmModule" )

INSTALL(TARGETS Integrator
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
