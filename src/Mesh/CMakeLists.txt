cmake_minimum_required(VERSION 2.8.5)

project(Mesh CXX)

include(charm)

add_library(Mesh
            DerivedData.C
            Reorder.C
            STLMesh.C
)

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
