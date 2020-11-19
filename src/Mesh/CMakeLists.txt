cmake_minimum_required(VERSION 2.8.5)

project(Mesh CXX)

include(charm)

if(HAS_EXAM2M)
  set(CHARMMESH CharmMesh.cpp)
endif()

add_library(Mesh
            DerivedData.cpp
            Gradients.cpp
            Reorder.cpp
            STLMesh.cpp
            ${CHARMMESH}
)

target_include_directories(Mesh PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PROJECT_BINARY_DIR}/../Mesh
                           ${PROJECT_BINARY_DIR}/../IO
                           ${TPL_INCLUDE_DIR}
                           ${PEGTL_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(Mesh PROPERTIES LIBRARY_OUTPUT_NAME quinoa_mesh)

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

if(HAS_EXAM2M)
  addCharmModule( "charmmesh" "Mesh" )
endif()
