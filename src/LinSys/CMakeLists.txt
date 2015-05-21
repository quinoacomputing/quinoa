cmake_minimum_required(VERSION 2.8.5)

project(LinSys CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(LinSys
            LinSysMerger.C
            ZoltanInterOp.C
)

addCharmModule( "linsysmerger" "LinSys" )
removeWarnings("ZoltanInterOp.C;LinSysMerger.C")

# Add extra dependency of LinSys on charm modules implementing advanced array
# element placement using maps. This is required to pull in the dependencies of
# the possible specializations of LinSys using various map placement strategies.
# For example, LinSys is specialized with Conductor which in turn depends on
# array placement maps.
add_dependencies( "LinSys" "linearmapCharmModule" )
add_dependencies( "LinSys" "unsmeshmapCharmModule" )

INSTALL(TARGETS LinSys
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
