## UnitTest executable ########################################################

if (HAS_MKL)
  set(TestMKLBetaMethod "Control/Options/TestMKLBetaMethod.C")
  set(TestMKLGammaMethod "Control/Options/TestMKLGammaMethod.C")
  set(TestMKLGaussianMethod "Control/Options/TestMKLGaussianMethod.C")
  set(TestMKLUniformMethod "Control/Options/TestMKLUniformMethod.C")
  set(TestMKLRNG "RNG/TestMKLRNG.C")
endif()

if(HAS_RNGSSE2)
  set(TestRNGSSE "RNG/TestRNGSSE.C")
endif()

add_executable(${UNITTEST_EXECUTABLE}
               UnitTestDriver.C
               UnitTest.C
               ../../tests/unit/Base/TestContainerUtil.C
               ../../tests/unit/Base/TestData.C
               ../../tests/unit/Base/TestException.C
               ../../tests/unit/Base/TestExceptionMPI.C
               ../../tests/unit/Base/TestFactory.C
               ../../tests/unit/Base/TestFlip_map.C
               ../../tests/unit/Base/TestHas.C
               ../../tests/unit/Base/TestPrint.C
               ../../tests/unit/Base/TestProcessControl.C
               ../../tests/unit/Base/TestPUPUtil.C
               ../../tests/unit/Base/TestReader.C
               ../../tests/unit/Base/TestStrConvUtil.C
               ../../tests/unit/Base/TestTaggedTuple.C
               ../../tests/unit/Base/TestTimer.C
               ../../tests/unit/Base/TestVector.C
               ../../tests/unit/Base/TestWriter.C
               ../../tests/unit/${TestMKLUniformMethod}
               ../../tests/unit/${TestMKLGaussianMethod}
               ../../tests/unit/${TestMKLBetaMethod}
               ../../tests/unit/${TestMKLGammaMethod}
               ../../tests/unit/Control/Options/TestRNG.C
               ../../tests/unit/Control/TestControl.C
               ../../tests/unit/Control/TestFileParser.C
               ../../tests/unit/Control/TestStringParser.C
               ../../tests/unit/Control/TestSystemComponents.C
               ../../tests/unit/Control/TestToggle.C
               ../../tests/unit/${TestScheme}
               ../../tests/unit/${TestError}
               ../../tests/unit/IO/TestExodusIIMeshReader.C
               ../../tests/unit/IO/TestMesh.C
               ../../tests/unit/IO/TestMeshReader.C
               ../../tests/unit/LoadBalance/TestLinearMap.C
               ../../tests/unit/LoadBalance/TestLoadDistributor.C
               ../../tests/unit/LoadBalance/TestUnsMeshMap.C
               ../../tests/unit/Mesh/TestAround.C
               ../../tests/unit/Mesh/TestDerivedData.C
               ../../tests/unit/Mesh/TestDerivedData_MPISingle.C
               ../../tests/unit/Mesh/TestGradients.C
               ../../tests/unit/Mesh/TestReorder.C
               ../../tests/unit/${TestMKLRNG}
               ../../tests/unit/${TestRNGSSE}
               ../../tests/unit/RNG/TestRNG.C
               ../../tests/unit/RNG/TestRandom123.C)

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
