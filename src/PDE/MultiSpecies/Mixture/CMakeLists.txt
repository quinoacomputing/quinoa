add_library(MultiSpeciesMixture
            Mixture.cpp)

target_include_directories(MultiSpeciesMixture PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/PDE
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main
                           ${CHARM_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${CBLAS_INCLUDE_DIRS})

set_target_properties(MultiSpeciesMixture
                      PROPERTIES LIBRARY_OUTPUT_NAME quinoa_multispeciesmixture)

INSTALL(TARGETS MultiSpeciesMixture
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
