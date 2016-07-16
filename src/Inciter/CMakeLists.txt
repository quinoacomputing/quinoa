cmake_minimum_required(VERSION 2.8.5)

project(Inciter CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Inciter
            Conductor.C
            Performer.C
            Tracker.C
            Partitioner.C
)

addCharmModule( "conductor" "Inciter" )
addCharmModule( "performer" "Inciter" )
addCharmModule( "tracker" "Inciter" )
addCharmModule( "partitioner" "Inciter" )

# Add extra dependency of Inciter on inciterCharmModule. This is required as
# Inciter refers to the main Charm++ proxy defined in the Charm++ module
# inciter (in Main/Inciter.C).
add_dependencies( "Inciter" "inciterCharmModule" )

## Add extra dependency of Conductor charm module on charm modules implementing
## advanced array element placement using maps. This is required since Conductor
## charm module depends on and uses advanced array element placement maps.
#add_dependencies( "conductorCharmModule" "linearmapCharmModule" )
#add_dependencies( "conductorCharmModule" "unsmeshmapCharmModule" )

# Add extra dependency of Inciter on LinSysMerger charm module. This is required
# so that linsysmerger.decl.h and linsysmerger.def.h are generated before
# Inciter including those.
add_dependencies( "Inciter" "linsysmergerCharmModule" )

INSTALL(TARGETS Inciter
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
