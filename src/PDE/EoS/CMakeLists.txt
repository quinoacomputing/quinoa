add_library(EOS
            EOS.cpp
            StiffenedGas.cpp
            JWL.cpp
            SmallShearSolid.cpp
            GodunovRomenskiAluminum.cpp)

target_include_directories(EOS PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/PDE
                           ${PROJECT_BINARY_DIR}/../Main
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${CBLAS_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(EOS
                      PROPERTIES LIBRARY_OUTPUT_NAME quinoa_eos)

INSTALL(TARGETS EOS
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
