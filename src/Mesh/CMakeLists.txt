cmake_minimum_required(VERSION 2.8.5)

project(Mesh CXX)

add_library(Mesh
            UnsMesh.C
            STLMesh.C
)

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
