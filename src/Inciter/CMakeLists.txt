cmake_minimum_required(VERSION 2.8.5)

project(Inciter CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Inciter
            Transporter.C
            Partitioner.C
            Carrier.C
            FluxCorrector.C
)

addCharmModule( "transporter" "Inciter" )
addCharmModule( "partitioner" "Inciter" )
addCharmModule( "carrier" "Inciter" )

# Add extra dependency of Inciter on inciterCharmModule. This is required as
# Inciter refers to the main Charm++ proxy defined in the Charm++ module
# inciter (in Main/Inciter.C).
add_dependencies( "Inciter" "inciterCharmModule" )

## Add extra dependency of Transporter charm module on charm modules implementing
## advanced array element placement using maps. This is required since
## Transporter charm module depends on and uses advanced array element placement
## maps.
#add_dependencies( "transporterCharmModule" "linearmapCharmModule" )
#add_dependencies( "transporterCharmModule" "unsmeshmapCharmModule" )

# Add extra dependency of Inciter on LinSysMerger charm module. This is required
# so that linsysmerger.decl.h and linsysmerger.def.h are generated before
# Inciter including those.
add_dependencies( "Inciter" "linsysmergerCharmModule" )

# Add extra dependency of Inciter on ParticleWriter charm module. This is
# required so that particlewriter.decl.h and particlewriter.def.h are generated
# before Inciter including those.
add_dependencies( "Inciter" "particlewriterCharmModule" )

# Add extra dependency of Carrier charm module on ParticleWriter charm module.
# This is required so that particlewriter.decl.h and particlewriter.def.h are
# generated before Carrier including those.
add_dependencies( "carrierCharmModule" "particlewriterCharmModule" )

# Add extra dependency of Transporter charm module on Partitioner charm module.
# This is required so that Partitioner.decl.h and Partitioner.def.h are
# generated before Transporter including those.
add_dependencies( "transporterCharmModule" "partitionerCharmModule" )

set_target_properties(Inciter PROPERTIES LIBRARY_OUTPUT_NAME quinoa_inciter)

INSTALL(TARGETS Inciter
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
