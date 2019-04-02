add_library(Beta
            BetaCoeffPolicy.C
            MassFractionBetaCoeffPolicy.C
            MixMassFractionBetaCoeffPolicy.C
            MixNumberFractionBetaCoeffPolicy.C
            NumberFractionBetaCoeffPolicy.C
            ConfigureBeta.C
            ConfigureNumberFractionBeta.C
            ConfigureMassFractionBeta.C
            ConfigureMixNumberFractionBeta.C
            ConfigureMixMassFractionBeta.C
)

target_include_directories(Beta PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/RNG
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Statistics
                           ${QUINOA_SOURCE_DIR}/DiffEq
                           ${QUINOA_SOURCE_DIR}/DiffEq/Velocity
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${MKL_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS})

set_target_properties(Beta PROPERTIES LIBRARY_OUTPUT_NAME quinoa_beta)

INSTALL(TARGETS Beta
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
