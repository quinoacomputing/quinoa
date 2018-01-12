cmake_minimum_required(VERSION 2.8.5)

project(Mesh CXX)

include(charm)

add_library(Mesh
            DerivedData.C
            Gradients.C
            Reorder.C
            STLMesh.C
)

set_target_properties(Mesh PROPERTIES LIBRARY_OUTPUT_NAME quinoa_mesh)

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
