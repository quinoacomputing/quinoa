cmake_minimum_required(VERSION 2.8.5)

project(Inciter CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Inciter
            Conductor.C
            Performer.C
            InciterSetup.C
)

addCharmModule( "conductor" "Inciter" )
addCharmModule( "performer" "Inciter" )
removeWarnings( "Conductor.C;Performer.C;InciterSetup.C" )

# Add extra dependency of Inciter on inciterCharmModule. This is required as
# Inciter refers to the main Charm++ proxy defined in the Charm++ module
# inciter (in Main/Inciter.C).
add_dependencies( "Inciter" "inciterCharmModule" )

INSTALL(TARGETS Inciter
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
