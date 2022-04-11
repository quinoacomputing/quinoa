add_library(MultiMatProblem
            UserDefined.cpp
            BoxInitialization.cpp
            FieldOutput.cpp
            SodShocktube.cpp
            InterfaceAdvection.cpp
            WaterAirShocktube.cpp
            ShockHeBubble.cpp
            UnderwaterEx.cpp
            ShockDensityWave.cpp
            EquilInterfaceAdvect.cpp)

target_include_directories(MultiMatProblem PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Control/Inciter
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/PDE
                           ${PROJECT_BINARY_DIR}/../Main
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(MultiMatProblem
                      PROPERTIES LIBRARY_OUTPUT_NAME quinoa_multimatproblem)

INSTALL(TARGETS MultiMatProblem
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
