cmake_minimum_required(VERSION 2.8.5)

project(Mesh CXX)

include(charm)

add_library(Mesh
            DerivedData.C
            Reorder.C
            STLMesh.C
)

removeWarnings( "DerivedData.C" )

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
