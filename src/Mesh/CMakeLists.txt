cmake_minimum_required(VERSION 2.8.5)

project(Mesh CXX)

include(charm)

add_library(Mesh
            DerivedData.C
            Gradients.C
            Reorder.C
            STLMesh.C
)

target_include_directories(Mesh PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${TPL_INCLUDE_DIR}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(Mesh PROPERTIES LIBRARY_OUTPUT_NAME quinoa_mesh)

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
