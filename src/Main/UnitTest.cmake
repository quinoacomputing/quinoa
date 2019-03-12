## UnitTest executable ########################################################

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
               ../UnitTest/tests/Mesh/TestDerivedData_MPISingle.C
               ../UnitTest/tests/Mesh/TestGradients.C
               ../UnitTest/tests/Mesh/TestReorder.C
               ../UnitTest/${TestMKLRNG}
               ../UnitTest/${TestRNGSSE}
               ../UnitTest/tests/RNG/TestRNG.C
               ../UnitTest/tests/RNG/TestRandom123.C)

target_include_directories(${UNITTEST_EXECUTABLE} PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/UnitTest
                           ${QUINOA_SOURCE_DIR}/LoadBalance
                           ${QUINOA_SOURCE_DIR}/IO
                           ${QUINOA_SOURCE_DIR}/RNG
                           ${TUT_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../UnitTest
                           ${PROJECT_BINARY_DIR}/../IO)

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
                      ${MKL_INTERFACE_LIBRARY}
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
