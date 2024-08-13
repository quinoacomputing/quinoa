project(Control CXX)

#### Quinoa control build configuration ########################################

include(charm)

#### Inciter control ###########################################################
if (ENABLE_INCITER)
  add_library(InciterControl
              StringParser.cpp
              Inciter/InputDeck/LuaParser.cpp
              Inciter/CmdLine/Parser.cpp)

  target_include_directories(InciterControl PUBLIC
                             ${QUINOA_SOURCE_DIR}
                             ${QUINOA_SOURCE_DIR}/Base
                             ${QUINOA_SOURCE_DIR}/Control
                             ${QUINOA_SOURCE_DIR}/Mesh
                             ${QUINOA_SOURCE_DIR}/PDE
                             ${QUINOA_SOURCE_DIR}/Inciter
                             ${PEGTL_INCLUDE_DIRS}
                             ${CHARM_INCLUDE_DIRS}
                             ${BRIGAND_INCLUDE_DIRS}
                             ${LUA_INCLUDE_DIR}
                             ${SOL2_INCLUDE_DIRS}
                             ${PROJECT_BINARY_DIR}/../Main)

  set_target_properties(InciterControl PROPERTIES
                        LIBRARY_OUTPUT_NAME quinoa_incitercontrol)

  INSTALL(TARGETS InciterControl
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
endif()

#### UnitTest control ##########################################################
if (ENABLE_UNITTEST)
  project(UnitTestControl CXX)

  add_library(UnitTestControl
              StringParser.cpp
              UnitTest/CmdLine/Parser.cpp)

  target_include_directories(UnitTestControl PUBLIC
                             ${QUINOA_SOURCE_DIR}
                             ${QUINOA_SOURCE_DIR}/Base
                             ${QUINOA_SOURCE_DIR}/Control
                             ${PROJECT_BINARY_DIR}/../Main
                             ${PEGTL_INCLUDE_DIRS}
                             ${CHARM_INCLUDE_DIRS}
                             ${BRIGAND_INCLUDE_DIRS}
                             ${OMEGA_H_INCLUDE_DIRS})

  set_target_properties(UnitTestControl PROPERTIES
                        LIBRARY_OUTPUT_NAME quinoa_unittestcontrol)

  INSTALL(TARGETS UnitTestControl
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
endif()

#### MeshConv control ##########################################################
if (ENABLE_MESHCONV)
  project(MeshConvControl CXX)

  add_library(MeshConvControl
              StringParser.cpp
              MeshConv/CmdLine/Parser.cpp)

  target_include_directories(MeshConvControl PUBLIC
                             ${QUINOA_SOURCE_DIR}
                             ${QUINOA_SOURCE_DIR}/Base
                             ${QUINOA_SOURCE_DIR}/IO
                             ${QUINOA_SOURCE_DIR}/Control
                             ${PROJECT_BINARY_DIR}/../Main
                             ${PEGTL_INCLUDE_DIRS}
                             ${CHARM_INCLUDE_DIRS}
                             ${BRIGAND_INCLUDE_DIRS})

  set_target_properties(MeshConvControl PROPERTIES
                        LIBRARY_OUTPUT_NAME quinoa_meshconvcontrol)

  INSTALL(TARGETS MeshConvControl
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
endif()
