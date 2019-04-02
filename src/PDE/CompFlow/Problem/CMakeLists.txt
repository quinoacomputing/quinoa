add_library(CompFlowProblem
            NLEnergyGrowth.C
            RayleighTaylor.C
            SodShocktube.C
            RotatedSodShocktube.C
            SedovBlastwave.C
            TaylorGreen.C
            UserDefined.C
            VorticalFlow.C)

target_include_directories(CompFlowProblem PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/PDE
                           ${PROJECT_BINARY_DIR}/../Main
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(CompFlowProblem
                      PROPERTIES LIBRARY_OUTPUT_NAME quinoa_compflowproblem)

INSTALL(TARGETS CompFlowProblem
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
