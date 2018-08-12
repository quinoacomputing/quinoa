cmake_minimum_required(VERSION 2.8.5)

project(DiffEq CXX)

add_library(DiffEq
            DiffEqStack.C
            ConfigureDirichlet.C
            ConfigureGeneralizedDirichlet.C
            ConfigureWrightFisher.C
            ConfigureOrnsteinUhlenbeck.C
            ConfigureDiagOrnsteinUhlenbeck.C
            ConfigureBeta.C
            ConfigureNumberFractionBeta.C
            ConfigureMassFractionBeta.C
            ConfigureMixNumberFractionBeta.C
            ConfigureMixMassFractionBeta.C
            ConfigureGamma.C
            ConfigureSkewNormal.C
            ConfigureVelocity.C
            ConfigurePosition.C
            ConfigureDissipation.C
)

set_target_properties(DiffEq PROPERTIES LIBRARY_OUTPUT_NAME quinoa_diffeq)

INSTALL(TARGETS DiffEq
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
