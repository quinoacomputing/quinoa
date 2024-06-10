project(quinoa CXX)

### Quinoa executables build configuration #####################################

# Basic configuration
add_definitions()

# Include function for adding Charm++ modules
include(charm)

# Link executables with the charmc wrapper
STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>"
       "${LINKER_COMPILER} -module CommonLBs -module collidecharm ${EXTRA_LINK_ARGS} -c++ <CMAKE_CXX_COMPILER>"
       CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")

include(ConfigExecutable)

# Configure configuration library
if (BRIGAND_FOUND)
  add_library(Config ${PROJECT_BINARY_DIR}/QuinoaConfig.cpp)

  set_target_properties(Config PROPERTIES LIBRARY_OUTPUT_NAME quinoa_config)

  INSTALL(TARGETS Config
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

  # Configure executable initialization library
  add_library(Init Init.cpp)

  set_target_properties(Init PROPERTIES LIBRARY_OUTPUT_NAME quinoa_init)

  target_include_directories(Init PUBLIC
                             ${QUINOA_SOURCE_DIR}
                             ${QUINOA_SOURCE_DIR}/Base
                             ${QUINOA_SOURCE_DIR}/Control
                             ${QUINOA_SOURCE_DIR}/Transfer
                             ${QUINOA_TPL_DIR}
                             ${PEGTL_INCLUDE_DIRS}
                             ${CHARM_INCLUDE_DIRS}
                             ${BRIGAND_INCLUDE_DIRS}
                             ${PROJECT_BINARY_DIR}/../Main)

  add_dependencies( "Init" "charestatecollectorCharmModule" )

  INSTALL(TARGETS Init
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
endif()

# Configure cmake variables for optional components in executables

if (ENABLE_INCITER)
  set(TestError "../../tests/unit/Inciter/AMR/TestError.cpp")
  set(TestScheme "../../tests/unit/Inciter/TestScheme.cpp")
  set(MESHREFINEMENT "MeshRefinement")
endif()

# Configure executable targets

if (ENABLE_UNITTEST AND ENABLE_TESTS)
  include("UnitTest.cmake")
endif()

if (ENABLE_INCITER)
  include("Inciter.cmake")
endif()

if (ENABLE_MESHCONV)
  include("MeshConv.cmake")
endif()
