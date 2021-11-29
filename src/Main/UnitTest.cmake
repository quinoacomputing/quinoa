## UnitTest executable ########################################################

add_executable(${UNITTEST_EXECUTABLE}
               UnitTestDriver.cpp
               UnitTest.cpp
               ../../tests/unit/Base/TestContainerUtil.cpp
               ../../tests/unit/Base/TestData.cpp
               ../../tests/unit/Base/TestException.cpp
               ../../tests/unit/Base/TestExceptionMPI.cpp
               ../../tests/unit/Base/TestFactory.cpp
               ../../tests/unit/Base/TestFlip_map.cpp
               ../../tests/unit/Base/TestHas.cpp
               ../../tests/unit/Base/TestPrint.cpp
               ../../tests/unit/Base/TestPUPUtil.cpp
               ../../tests/unit/Base/TestReader.cpp
               ../../tests/unit/Base/TestPrintUtil.cpp
               ../../tests/unit/Base/TestTaggedTuple.cpp
               ../../tests/unit/Base/TestTaggedTuplePrint.cpp
               ../../tests/unit/Base/TestTaggedTupleDeepPrint.cpp
               ../../tests/unit/Base/TestTimer.cpp
               ../../tests/unit/Base/TestVector.cpp
               ../../tests/unit/Base/TestWriter.cpp
               ../../tests/unit/Control/TestFileParser.cpp
               ../../tests/unit/Control/TestStringParser.cpp
               ../../tests/unit/Control/TestSystemComponents.cpp
               ../../tests/unit/Control/TestToggle.cpp
               ../../tests/unit/${TestScheme}
               ../../tests/unit/${TestError}
               ../../tests/unit/IO/TestExodusIIMeshReader.cpp
               ../../tests/unit/IO/TestMesh.cpp
               ../../tests/unit/IO/TestMeshReader.cpp
               ../../tests/unit/LoadBalance/TestLinearMap.cpp
               ../../tests/unit/LoadBalance/TestLoadDistributor.cpp
               ../../tests/unit/LoadBalance/TestUnsMeshMap.cpp
               ../../tests/unit/LinearSolver/TestCSR.cpp
               ../../tests/unit/LinearSolver/TestConjugateGradients.cpp
               ../../tests/unit/Mesh/TestAround.cpp
               ../../tests/unit/Mesh/TestDerivedData.cpp
               ../../tests/unit/Mesh/TestDerivedData_MPISingle.cpp
               ../../tests/unit/Mesh/TestGradients.cpp
               ../../tests/unit/Mesh/TestReorder.cpp)

target_include_directories(${UNITTEST_EXECUTABLE} PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/UnitTest
                           ${QUINOA_SOURCE_DIR}/LoadBalance
                           ${QUINOA_SOURCE_DIR}/IO
                           ${TUT_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${RANDOM123_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../UnitTest
                           ${PROJECT_BINARY_DIR}/../IO)

config_executable(${UNITTEST_EXECUTABLE})

target_link_libraries(${UNITTEST_EXECUTABLE}
                      Base
                      Config
                      Init
                      ${MESHREFINEMENT}
                      UnitTest
                      UnitTestControl
                      LoadBalance
                      LinearSolver
                      Mesh
                      MeshDetect
                      NativeMeshIO
                      ExodusIIMeshIO
                      HyperMeshIO
                      ${PUGIXML_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${HDF5_HL_LIBRARIES}      # only for static link
                      ${HDF5_C_LIBRARIES}
                      ${BACKWARD_LIBRARIES}
                      ${LIBCXX_LIBRARIES}       # only for static link with libc++
                      ${LIBCXXABI_LIBRARIES})   # only for static link with libc++

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
