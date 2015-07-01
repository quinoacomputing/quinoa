

SET( unit_tests
     ArrayMeshTest 
     AspectRatioGammaTest 
     AveragingQMTest
     AWMetricTest
     AWQualityMetricTest
     BCDTest
     BoundedCylinderDomainTest
     CachingTargetTest
     CircleDomainTest
     CompositeMetricTest
     CompositeOFTest
     ConditionNumberTest
     ConicDomainTest
     CylinderDomainTest
     DomainClassifierTest
     EdgeLengthMetricTest
     EdgeLengthQualityMetricTest
     ElemSampleQMTest
     ExtraDataTest
     FileTokenizerTest
     GeomPrimTest
     HexLagrangeShapeTest
     HigherOrderTest
     IdealElementTest
     IdealTargetTest
     IdealWeightMeanRatioTest
     IdealWeightInverseMeanRatioTest
     InstructionQueueTest
     LinearMappingFunctionTest
     LineDomainTest
     LVQDTargetTest
     MappingFunctionTest
     Matrix3DTest
     MeshInterfaceTest
     MeshImplTest
     MiscTargetMetricTests
     MsqFreeVertexIndexIteratorTest
     MsqHessianTest
     MsqMatrixTest
     MsqMeshEntityTest
     MsqVertexTest
     NodeSetTest
     NumericalOFTest
     ObjectiveFunctionTest
     PatchDataTest
     PatchDataTestNormals
     PatchSetTest
     PlanarGeometryTest
     PMeanPMetricTest
     PMeanPTemplateTest
     QuadLagrangeShapeTest
     QualityAssessorTest
     QualityMetricTest
     ShapeOrientTargetMetricTests
     ShapeSizeOrientTargetMetricTests
     ShapeSizeTargetMetricTests
     ShapeTargetMetricTests
     SizeTargetMetricTests
     SlaveBoundaryVerticesTest
     SphericalDomainTest
     SphericalGeometryTest
     StdDevTemplateTest
     SymMatrix3DTest
     TagVertexMeshTest
     TargetCalculatorTest
     TargetReadWriteTest
     TerminationCriterionTest
     TetLagrangeShapeTest
     TriLagrangeShapeTest
     TopologyInfoTest
     TMPDerivsTest
     TMetricTest
     TQualityMetricTest
     UntangleBetaTest
     UntangleTargetMetricTests
     Vector3DTest
     VertexConditionNumberTest
     VertexCullingRegressionTest
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
