add_library(Dirichlet
            DirichletCoeffPolicy.cpp
            GeneralizedDirichletCoeffPolicy.cpp
            MixDirichletCoeffPolicy.cpp
            ConfigureDirichlet.cpp
            ConfigureMixDirichlet.cpp
            ConfigureGeneralizedDirichlet.cpp
)

target_include_directories(Dirichlet PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/RNG
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Statistics
                           ${QUINOA_SOURCE_DIR}/DiffEq
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${MKL_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS})

set_target_properties(Dirichlet PROPERTIES LIBRARY_OUTPUT_NAME quinoa_dirichlet)

INSTALL(TARGETS Dirichlet
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
