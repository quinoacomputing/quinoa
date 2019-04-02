add_library(Gamma
            GammaCoeffPolicy.C
            ConfigureGamma.C
)

target_include_directories(Gamma PUBLIC
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

set_target_properties(Gamma PROPERTIES LIBRARY_OUTPUT_NAME quinoa_gamma)

INSTALL(TARGETS Gamma
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
