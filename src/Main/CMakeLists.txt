cmake_minimum_required(VERSION 2.8.5)

project(quinoa CXX)

### Quinoa executables build configuration #####################################

# Basic configuration
add_definitions()

# Include function for adding Charm++ modules
include(charm)

# Link executables with the charmc wrapper
STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>"
       ${CHARM_COMPILER}
       CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")

### Inciter executable #########################################################
add_executable(${INCITER_EXECUTABLE}
               InciterDriver.C
               InciterPrint.C
               Inciter.C
               #Inciter_main.C
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${INCITER_EXECUTABLE} PROPERTIES
                        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${INCITER_EXECUTABLE} PROPERTIES
                        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${INCITER_EXECUTABLE}
                      Base
                      IO
                      Mesh
                      Statistics
                      InciterControl
                      ${NEMESIS_LIBRARIES}
                      ${EXODUS_LIBRARIES}
                      ${NETCDF_LIBRARIES}
                      ${HDF5_LIBRARY}
                      ${ZOLTAN_LIBRARIES}
)

# Add custom dependencies for Inciter's main Charm++ module
addCharmModule( "inciter" "${INCITER_EXECUTABLE}" )
removeWarnings( "InciterDriver.C;InciterPrint.C;Inciter.C")

INSTALL(TARGETS ${INCITER_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

### RNGTest executable #########################################################
add_executable(${RNGTEST_EXECUTABLE}
               RNGTestDriver.C
               RNGPrint.C
               RNGTest.C
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${RNGTEST_EXECUTABLE} PROPERTIES
                        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${RNGTEST_EXECUTABLE} PROPERTIES
                        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${RNGTEST_EXECUTABLE}
                      Base
                      RNG
                      RNGTest
                      RNGTestControl
                      ${TESTU01_PROBDIST_LIBRARY}
                      ${TESTU01_LIBRARY}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${RNGSSE_LIBRARY}
)

# Add custom dependencies for RNGTest's main Charm++ module
addCharmModule( "rngtest" "${RNGTEST_EXECUTABLE}" )
removeWarnings( "RNGTest.C;RNGTestDriver.C;RNGPrint.C")

INSTALL(TARGETS ${RNGTEST_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

## UnitTest executable ########################################################
add_executable(${UNITTEST_EXECUTABLE}
               UnitTestDriver.C
               UnitTest.C
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${UNITTEST_EXECUTABLE} PROPERTIES
                        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${UNITTEST_EXECUTABLE} PROPERTIES
                        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${UNITTEST_EXECUTABLE}
                      Base
                      IO
                      Mesh
                      UnitTest
                      UnitTestControl
                      #${SILO_LIBRARY}
                      ${NEMESIS_LIBRARIES}
                      ${EXODUS_LIBRARIES}
                      ${NETCDF_LIBRARIES}
                      ${ZOLTAN_LIBRARIES}
)

# Add custom dependencies for UnitTest's main Charm++ module
addCharmModule( "unittest" "${UNITTEST_EXECUTABLE}" )
removeWarnings( "UnitTest.C;UnitTestDriver.C")

INSTALL(TARGETS ${UNITTEST_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

## RegTest executable #########################################################
add_executable(${REGTEST_EXECUTABLE}
               RegTestDriver.C
               RegTest.C
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${REGTEST_EXECUTABLE} PROPERTIES
                        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${REGTEST_EXECUTABLE} PROPERTIES
                        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${REGTEST_EXECUTABLE}
                      Base
                      RegTest
                      RegTestControl
)

# Add custom dependencies for RegTest's main Charm++ module
addCharmModule( "regtest" "${REGTEST_EXECUTABLE}" )
removeWarnings( "RegTest.C;RegTestDriver.C")

INSTALL(TARGETS ${REGTEST_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

### MeshConv executable ########################################################
add_executable(${MESHCONV_EXECUTABLE}
               MeshConvDriver.C
               MeshConv.C
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${MESHCONV_EXECUTABLE} PROPERTIES
                        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${MESHCONV_EXECUTABLE} PROPERTIES
                        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${MESHCONV_EXECUTABLE}
                      Base
                      IO
                      Mesh
                      MeshConvControl
                      #${SILO_LIBRARY}
                      ${NEMESIS_LIBRARIES}
                      ${EXODUS_LIBRARIES}
                      ${NETCDF_LIBRARIES}
                      ${ZOLTAN_LIBRARIES}
)

# Add custom dependencies for MeshConv's main Charm++ module
addCharmModule( "meshconv" "${MESHCONV_EXECUTABLE}" )
removeWarnings( "MeshConv.C;MeshConvDriver.C" )

INSTALL(TARGETS ${MESHCONV_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

### Walker executable ##########################################################
add_executable(${WALKER_EXECUTABLE}
               WalkerDriver.C
               WalkerPrint.C
               Walker.C
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${WALKER_EXECUTABLE} PROPERTIES
                        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${WALKER_EXECUTABLE} PROPERTIES
                        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${WALKER_EXECUTABLE}
                      Base
                      DiffEq
                      IO
                      Mesh
                      Statistics
                      RNG
                      Walker
                      WalkerControl
                      ${NEMESIS_LIBRARIES}
                      ${EXODUS_LIBRARIES}
                      ${NETCDF_LIBRARIES}
                      ${HDF5_LIBRARY}
                      ${LAPACKE_LIB}
                      ${LAPACK_LIBRARIES}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${RNGSSE_LIBRARY}
                      ${ZOLTAN_LIBRARIES}
)

# Add custom dependencies for Walker's main Charm++ module
addCharmModule( "walker" "${WALKER_EXECUTABLE}" )
removeWarnings( "WalkerDriver.C;WalkerPrint.C;Walker.C")

INSTALL(TARGETS ${WALKER_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

