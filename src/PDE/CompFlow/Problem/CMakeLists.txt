add_library(CompFlowProblem
            FieldOutput.cpp
            NLEnergyGrowth.cpp
            RayleighTaylor.cpp
            SodShocktube.cpp
            SheddingFlow.cpp
            RotatedSodShocktube.cpp
            SedovBlastwave.cpp
            TaylorGreen.cpp
            UserDefined.cpp
            VorticalFlow.cpp
            GaussHumpCompflow.cpp
            ShockDensityWave.cpp)

target_include_directories(CompFlowProblem PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Control/Inciter
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/PDE
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${PROJECT_BINARY_DIR}/../Main
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${CBLAS_INCLUDE_DIRS})

set_target_properties(CompFlowProblem
                      PROPERTIES LIBRARY_OUTPUT_NAME quinoa_compflowproblem)

INSTALL(TARGETS CompFlowProblem
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
