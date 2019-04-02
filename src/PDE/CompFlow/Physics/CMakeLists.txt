add_library(CGCompFlowPhysics
            CGNavierStokes.C)

target_include_directories(CGCompFlowPhysics PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${PROJECT_BINARY_DIR}/../Main
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(CGCompFlowPhysics
                      PROPERTIES LIBRARY_OUTPUT_NAME quinoa_cgcompflowphysics)

INSTALL(TARGETS CGCompFlowPhysics
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
