add_library(Integrate
            Quadrature.C
            Initialize.C
            Mass.C
            Surface.C
            Boundary.C
            Volume.C
            Source.C
            Basis.C
            Riemann/RiemannFactory.C)

target_include_directories(Integrate PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/PDE
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${PROJECT_BINARY_DIR}/../Main
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(Integrate
                      PROPERTIES LIBRARY_OUTPUT_NAME quinoa_integrate)

INSTALL(TARGETS Integrate
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
