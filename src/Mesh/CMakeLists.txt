project(Mesh CXX)

add_library(Mesh
            DerivedData.cpp
            Gradients.cpp
            Reorder.cpp
            CommMap.cpp
            STLMesh.cpp)

target_include_directories(Mesh PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${PROJECT_BINARY_DIR}/../Main
                           ${TPL_INCLUDE_DIR}
                           ${PEGTL_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${CBLAS_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(Mesh PROPERTIES LIBRARY_OUTPUT_NAME quinoa_mesh)

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
