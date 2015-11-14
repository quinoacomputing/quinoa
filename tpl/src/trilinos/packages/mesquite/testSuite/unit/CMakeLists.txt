
# Tests that are in the following list depend on the nonexistant header file
# <Mesquite_cppunit/extensions/HelperMacros.h>.
# TODO ressurect this file, reactivate the tests
SET(unit_tests_helpermacros
     AspectRatioGammaTest
     AveragingQMTest
     BCDTest
     CachingTargetTest
     CompositeMetricTest
     CompositeOFTest
     ConditionNumberTest
     EdgeLengthMetricTest
     EdgeLengthQualityMetricTest
     FileTokenizerTest
     HexLagrangeShapeTest
     HigherOrderTest
     IdealElementTest
     IdealTargetTest
     IdealWeightMeanRatioTest
     IdealWeightInverseMeanRatioTest
     Matrix3DTest
     MeshInterfaceTest
     MsqFreeVertexIndexIteratorTest
     MsqHessianTest
     MsqMatrixTest
     MsqMeshEntityTest
     MsqVertexTest
     NumericalOFTest
     ObjectiveFunctionTest
     PatchDataTest
     PlanarGeometryTest
     PMeanPMetricTest
     QualityMetricTest
     SphericalGeometryTest
     TargetReadWriteTest
     TerminationCriterionTest
     TopologyInfoTest
     UntangleBetaTest
     Vector3DTest
     VertexConditionNumberTest
     VertexCullingRegressionTest
     )

SET( unit_tests
     ArrayMeshTest
     AWMetricTest
     AWQualityMetricTest
     BoundedCylinderDomainTest
     CircleDomainTest
     ConicDomainTest
     CylinderDomainTest
     DomainClassifierTest
     ElemSampleQMTest
     ExtraDataTest
     GeomPrimTest
     InstructionQueueTest
     LinearMappingFunctionTest
     LineDomainTest
     LVQDTargetTest
     MappingFunctionTest
     MeshImplTest
     MiscTargetMetricTests
     NodeSetTest
     PatchDataTestNormals
     PatchSetTest
     PMeanPTemplateTest
     QuadLagrangeShapeTest
     QualityAssessorTest
     ShapeOrientTargetMetricTests
     ShapeSizeOrientTargetMetricTests
     ShapeSizeTargetMetricTests
     ShapeTargetMetricTests
     SizeTargetMetricTests
     SlaveBoundaryVerticesTest
     SphericalDomainTest
     StdDevTemplateTest
     SymMatrix3DTest
     TagVertexMeshTest
     TargetCalculatorTest
     TetLagrangeShapeTest
     TriLagrangeShapeTest
     TMPDerivsTest
     TMetricTest
     TQualityMetricTest
     UntangleTargetMetricTests
     VtkTest
     XYRectangleTest )

SET( unit_libs ${CPPUNIT_LIBRARY} )

IF (${PACKAGE_NAME}_ENABLE_IMESH)
  SET(unit_tests ${unit_tests} iMeshTest)
  SET(unit_libs ${unit_libs} ${IMESH_LIBRARY})
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_ExodusII)
  SET(unit_tests ${unit_tests} ExodusTest)
  SET(unit_libs ${unit_libs} ${EXODUS_LIBRARY})
ENDIF()

SET( unit_srcs
     msq_test_main.cpp 
     QualityMetricTester.cpp 
     MesquiteTestRunner.hpp 
     ObjectiveFunctionTests.hpp 
     ObjectiveFunctionTests.cpp 
     PatchDataInstances.hpp 
     QualityMetricTester.hpp 
     UnitUtil.hpp 
     TargetMetricTest.hpp
     TMPQualityMetricTest.hpp )

FOREACH( test ${unit_tests} )
  SET( unit_srcs ${unit_srcs} ${test}.cpp )
ENDFOREACH()

IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")
  TRIBITS_ADD_EXECUTABLE( msq_test SOURCES ${unit_srcs} DEPLIBS ${CPPUNIT_LIBRARY} )
#  TARGET_LINK_LIBRARIES( msq_test mesquite ${CPPUNIT_LIBRARY} )
  FOREACH( test ${unit_tests} )
    TRIBITS_ADD_TEST( msq_test NAME ${test} ARGS ${test} NUM_MPI_PROCS 1 )
  ENDFOREACH()
ELSE()
  IF (${PACKAGE_NAME}_ENABLE_TESTS)
    ADD_EXECUTABLE( msq_test ${unit_srcs} )
    TARGET_LINK_LIBRARIES( msq_test mesquite ${CPPUNIT_LIBRARY} )
    FOREACH( test ${unit_tests} )
      ADD_TEST( ${test} ${EXECUTABLE_OUTPUT_PATH}/msq_test ${test} )
    ENDFOREACH()
  ENDIF()
ENDIF()

INCLUDE_DIRECTORIES( ${CPPUNIT_INCLUDES} )
SET_SOURCE_FILES_PROPERTIES(${unit_srcs} COMPILE_FLAGS "${TEST_DEFINES} -DSRCDIR=\"\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"\"")
