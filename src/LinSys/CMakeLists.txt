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

INSTALL(TARGETS LinSys
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
