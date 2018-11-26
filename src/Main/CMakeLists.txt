cmake_minimum_required(VERSION 2.8.5)

project(quinoa CXX)

### Quinoa executables build configuration #####################################

# Basic configuration
add_definitions()

# Include function for adding Charm++ modules
include(charm)

# Link executables with the charmc wrapper
STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>"
       "${LINKER_COMPILER} -module CommonLBs ${EXTRA_LINK_ARGS} -c++ <CMAKE_CXX_COMPILER>"
       CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")

include(ConfigExecutable)

# Configure configuration library
add_library(Config ${PROJECT_BINARY_DIR}/QuinoaConfig.C)

set_target_properties(Config PROPERTIES LIBRARY_OUTPUT_NAME quinoa_config)

INSTALL(TARGETS Config
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

# Configure executable initialization library
add_library(Init Init.C)

set_target_properties(Init PROPERTIES LIBRARY_OUTPUT_NAME quinoa_init)

target_include_directories(Init PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
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

if (HAS_ROOT)
  set(ROOTMESHIO "RootMeshIO")
endif()

if (HAS_OMEGA_H)
  set(OMEGAHMESHIO "OmegaHMeshIO")
endif()

if (ENABLE_INCITER)
  set(TestError "../UnitTest/tests/Inciter/AMR/TestError.C")
  set(TestScheme "../UnitTest/tests/Inciter/TestScheme.C")
  set(MESHREFINEMENT "MeshRefinement")
endif()

## UnitTest executable ########################################################

if (ENABLE_UNITTEST)

  if (HAS_MKL)
    set(TestMKLBetaMethod "tests/Control/Options/TestMKLBetaMethod.C")
    set(TestMKLGammaMethod "tests/Control/Options/TestMKLGammaMethod.C")
    set(TestMKLGaussianMethod "tests/Control/Options/TestMKLGaussianMethod.C")
    set(TestMKLUniformMethod "tests/Control/Options/TestMKLUniformMethod.C")
    set(TestMKLRNG "tests/RNG/TestMKLRNG.C")
  endif()

  if(HAS_RNGSSE2)
    set(TestRNGSSE "tests/RNG/TestRNGSSE.C")
  endif()

  add_executable(${UNITTEST_EXECUTABLE}
                 UnitTestDriver.C
                 UnitTest.C
                 ../UnitTest/tests/Base/TestContainerUtil.C
                 ../UnitTest/tests/Base/TestData.C
                 ../UnitTest/tests/Base/TestException.C
                 ../UnitTest/tests/Base/TestExceptionMPI.C
                 ../UnitTest/tests/Base/TestFactory.C
                 ../UnitTest/tests/Base/TestFlip_map.C
                 ../UnitTest/tests/Base/TestHas.C
                 ../UnitTest/tests/Base/TestPrint.C
                 ../UnitTest/tests/Base/TestProcessControl.C
                 ../UnitTest/tests/Base/TestPUPUtil.C
                 ../UnitTest/tests/Base/TestReader.C
                 ../UnitTest/tests/Base/TestStrConvUtil.C
                 ../UnitTest/tests/Base/TestTaggedTuple.C
                 ../UnitTest/tests/Base/TestTimer.C
                 ../UnitTest/tests/Base/TestVector.C
                 ../UnitTest/tests/Base/TestWriter.C
                 ../UnitTest/${TestMKLUniformMethod}
                 ../UnitTest/${TestMKLGaussianMethod}
                 ../UnitTest/${TestMKLBetaMethod}
                 ../UnitTest/${TestMKLGammaMethod}
                 ../UnitTest/tests/Control/Options/TestRNG.C
                 ../UnitTest/tests/Control/TestControl.C
                 ../UnitTest/tests/Control/TestFileParser.C
                 ../UnitTest/tests/Control/TestStringParser.C
                 ../UnitTest/tests/Control/TestSystemComponents.C
                 ../UnitTest/tests/Control/TestToggle.C
                 ../UnitTest/${TestScheme}
                 ../UnitTest/${TestError}
                 ../UnitTest/tests/IO/TestExodusIIMeshReader.C
                 ../UnitTest/tests/IO/TestMesh.C
                 ../UnitTest/tests/IO/TestMeshReader.C
                 ../UnitTest/tests/LoadBalance/TestLinearMap.C
                 ../UnitTest/tests/LoadBalance/TestLoadDistributor.C
                 ../UnitTest/tests/LoadBalance/TestUnsMeshMap.C
                 ../UnitTest/tests/Mesh/TestAround.C
                 ../UnitTest/tests/Mesh/TestDerivedData.C
                 ../UnitTest/tests/Mesh/TestGradients.C
                 ../UnitTest/tests/Mesh/TestReorder.C
                 ../UnitTest/${TestMKLRNG}
                 ../UnitTest/${TestRNGSSE}
                 ../UnitTest/tests/RNG/TestRNG.C
                 ../UnitTest/tests/RNG/TestRandom123.C)

  target_include_directories(${UNITTEST_EXECUTABLE} PUBLIC
                             ${QUINOA_SOURCE_DIR}
                             ${QUINOA_SOURCE_DIR}/UnitTest
                             ${QUINOA_SOURCE_DIR}/LinSys
                             ${QUINOA_SOURCE_DIR}/LoadBalance
                             ${QUINOA_SOURCE_DIR}/IO
                             ${QUINOA_SOURCE_DIR}/RNG
                             ${TUT_INCLUDE_DIRS}
                             ${LAPACKE_INCLUDE_DIRS}
                             ${PROJECT_BINARY_DIR}/../LinSys
                             ${PROJECT_BINARY_DIR}/../UnitTest)

  config_executable(${UNITTEST_EXECUTABLE})

  target_link_libraries(${UNITTEST_EXECUTABLE}
                        Base
                        Config
                        Init
                        RNG
                        ${MESHREFINEMENT}
                        UnitTest
                        UnitTestControl
                        LoadBalance
                        Mesh
                        MeshDetect
                        NativeMeshIO
                        ExodusIIMeshIO
                        HyperMeshIO
                        ${OMEGAHMESHIO}
                        ${PUGIXML_LIBRARIES}
                        ${SEACASExodus_LIBRARIES}
                        ${RNGSSE2_LIBRARIES}
                        ${MKL_INTERFACE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${MKL_CORE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${NETCDF_LIBRARIES}       # only for static link
                        ${HDF5_HL_LIBRARIES}      # only for static link
                        ${HDF5_C_LIBRARIES}
                        ${AEC_LIBRARIES}          # only for static link
                        ${BACKWARD_LIBRARIES}
                        ${OMEGA_H_LIBRARIES})

  # Add custom dependencies for UnitTest's main Charm++ module
  if(ENABLE_INCITER)
    addCharmModule( "unittestinciter" "${UNITTEST_EXECUTABLE}" )
    add_dependencies( "unittestinciterCharmModule"
                      "charestatecollectorCharmModule"
                      "mpirunnerinciterCharmModule" )
  else()
    addCharmModule( "unittest" "${UNITTEST_EXECUTABLE}" )
    add_dependencies( "unittestCharmModule"
                      "charestatecollectorCharmModule"
                      "mpirunnerCharmModule" )
  endif()

endif(ENABLE_UNITTEST)

### Inciter executable #########################################################

if (ENABLE_INCITER)

  add_executable(${INCITER_EXECUTABLE}
                 InciterDriver.C
                 InciterPrint.C
                 Inciter.C)

  config_executable(${INCITER_EXECUTABLE})

  target_link_libraries(${INCITER_EXECUTABLE}
                        PDE
                        InciterControl
                        Inciter
                        MeshRefinement
                        LoadBalance
                        ZoltanInterOp
                        LinSys
                        Base
                        Config
                        Init
                        IO
                        ExodusIIMeshIO
                        ${OMEGAHMESHIO}
                        ${ROOTMESHIO}
                        MeshDetect
                        Mesh
                        Statistics
                        ${SEACASExodus_LIBRARIES}
                        ${HYPRE_LIBRARIES}
                        ${ROOT_LIBRARIES}
                        ${Zoltan2_LIBRARIES}
                        ${RNGSSE2_LIBRARIES}
                        ${LAPACKE_LIBRARIES}      # only if MKL not found
                        ${MKL_INTERFACE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${MKL_CORE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${NETCDF_LIBRARIES}       # only for static link
                        ${HDF5_HL_LIBRARIES}      # only for static link
                        ${HDF5_C_LIBRARIES}
                        ${AEC_LIBRARIES}          # only for static link
                        ${BACKWARD_LIBRARIES}
                        ${OMEGA_H_LIBRARIES})

  # Add custom dependencies for Inciter's main Charm++ module
  addCharmModule( "inciter" "${INCITER_EXECUTABLE}" )

  add_dependencies( "inciterCharmModule" "charestatecollectorCharmModule" )

endif(ENABLE_INCITER)

### RNGTest executable #########################################################
if (ENABLE_RNGTEST)

  add_executable(${RNGTEST_EXECUTABLE}
                 RNGTestDriver.C
                 RNGPrint.C
                 RNGTest.C)

  target_include_directories(${RNGTEST_EXECUTABLE} PUBLIC
                             ${QUINOA_SOURCE_DIR}/RNGTest
                             ${QUINOA_SOURCE_DIR}/Main)

  config_executable(${RNGTEST_EXECUTABLE})

  target_link_libraries(${RNGTEST_EXECUTABLE}
                        RNG
                        RNGTest
                        RNGTestControl
                        Base
                        Config
                        Init
                        ${TESTU01_LIBRARIES}
                        ${MKL_INTERFACE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${MKL_CORE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${RNGSSE2_LIBRARIES}
                        ${BACKWARD_LIBRARIES})

  # Add custom dependencies for RNGTest's main Charm++ module
  addCharmModule( "rngtest" "${RNGTEST_EXECUTABLE}" )

  add_dependencies( "rngtestCharmModule" "charestatecollectorCharmModule" )

endif(ENABLE_RNGTEST)

### MeshConv executable ########################################################
if (ENABLE_MESHCONV)

  add_executable(${MESHCONV_EXECUTABLE}
                 MeshConvDriver.C
                 MeshConv.C)

  config_executable(${MESHCONV_EXECUTABLE})

  target_include_directories(${MESHCONV_EXECUTABLE} PUBLIC
                             ${PROJECT_BINARY_DIR}/../Base)

  target_link_libraries(${MESHCONV_EXECUTABLE}
                        NativeMeshIO
                        ExodusIIMeshIO
                        HyperMeshIO
                        MeshDetect
                        Mesh
                        MeshConvControl
                        Base
                        Config
                        Init
                        ${PUGIXML_LIBRARIES}
                        ${SEACASExodus_LIBRARIES}
                        ${NETCDF_LIBRARIES}       # only for static link
                        ${HDF5_HL_LIBRARIES}      # only for static link
                        ${HDF5_C_LIBRARIES}
                        ${AEC_LIBRARIES}          # only for static link
                        ${BACKWARD_LIBRARIES}
                        ${OMEGA_H_LIBRARIES})

  # Add custom dependencies for MeshConv's main Charm++ module
  addCharmModule( "meshconv" "${MESHCONV_EXECUTABLE}" )

  add_dependencies( "meshconvCharmModule" "charestatecollectorCharmModule" )

endif(ENABLE_MESHCONV)

### Walker executable ##########################################################
if (ENABLE_WALKER)

  add_executable(${WALKER_EXECUTABLE}
                 WalkerDriver.C
                 WalkerPrint.C
                 Walker.C)

  config_executable(${WALKER_EXECUTABLE})

  target_link_libraries(${WALKER_EXECUTABLE}
                        DiffEq
                        RNG
                        Walker
                        Statistics
                        IO
                        WalkerControl
                        Base
                        Config
                        Init
                        ${SEACASExodus_LIBRARIES}
                        ${LAPACKE_LIBRARIES}      # only if MKL not found
                        ${MKL_INTERFACE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${MKL_CORE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${RNGSSE2_LIBRARIES}
                        ${NETCDF_LIBRARIES}       # only for static link
                        ${HDF5_HL_LIBRARIES}      # only for static link
                        ${HDF5_C_LIBRARIES}
                        ${AEC_LIBRARIES}          # only for static link
                        ${BACKWARD_LIBRARIES})

  # Add custom dependencies for Walker's main Charm++ module
  addCharmModule( "walker" "${WALKER_EXECUTABLE}" )

  add_dependencies( "walkerCharmModule" "charestatecollectorCharmModule" )

endif(ENABLE_WALKER)

### FileConv executable ########################################################
if (ENABLE_FILECONV)

  add_executable(${FILECONV_EXECUTABLE}
                 FileConvDriver.C
                 FileConv.C)

  target_include_directories(${FILECONV_EXECUTABLE} PUBLIC
                             ${QUINOA_SOURCE_DIR}/IO
                             ${PROJECT_BINARY_DIR}/../Base)

  config_executable(${FILECONV_EXECUTABLE})

  target_link_libraries(${FILECONV_EXECUTABLE}
                        ExodusIIMeshIO
                        ${ROOTMESHIO}
                        Mesh
                        FileConvControl
                        Base
                        Config
                        Init
                        ${SEACASExodus_LIBRARIES}
                        ${ROOT_LIBRARIES}
                        ${NETCDF_LIBRARIES}       # only for static link
                        ${BACKWARD_LIBRARIES}
                        ${OMEGA_H_LIBRARIES})

  # Add custom dependencies for FileConv's main Charm++ module
  addCharmModule( "fileconv" "${FILECONV_EXECUTABLE}" )

  add_dependencies( "fileconvCharmModule" "charestatecollectorCharmModule" )

endif(ENABLE_FILECONV)
