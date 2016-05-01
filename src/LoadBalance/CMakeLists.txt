cmake_minimum_required(VERSION 2.8.5)

project(LoadBalance CXX)

include(charm)

add_library(LoadBalance
            CommMap.C
            LinearMap.C
            UnsMeshMap.C
)

addCharmModule( "linearmap" "LoadBalance" )
addCharmModule( "unsmeshmap" "LoadBalance" )
removeWarnings( "CommMap.C;LinearMap.C;UnsMeshMap.C")

INSTALL(TARGETS LoadBalance
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
