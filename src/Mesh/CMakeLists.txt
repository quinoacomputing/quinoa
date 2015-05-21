cmake_minimum_required(VERSION 2.8.5)

project(Mesh CXX)

include(charm)

add_library(Mesh
            DerivedData.C
            Reorder.C
            CommMap.C
            STLMesh.C
            UnsMeshMap.C
)

addCharmModule( "unsmeshmap" "Mesh" )
removeWarnings( "DerivedData.C;CommMap.C;UnsMeshMap.C" )

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
