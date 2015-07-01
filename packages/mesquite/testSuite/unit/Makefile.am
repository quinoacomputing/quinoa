@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/src/libmesquite.la
AM_LDFLAGS = @AM_LDFLAGS@
AM_CPPFLAGS = @AM_CPPFLAGS@


if WITH_CPPUNIT
  check_PROGRAMS = msq_test
  msq_test_DEPENDENCIES = $(top_builddir)/src/libmesquite.la
  TESTS = $(check_PROGRAMS)
  CPPUNIT_CPPFLAGS_SET = @CPPUNIT_CPPFLAGS@
  CPPUNIT_LDFLAGS_SET = @CPPUNIT_LDFLAGS@
else
  TESTS = no_cpp_unit
  msq_test_DEPENDENCIES = 
  CPPUNIT_CPPFLAGS_SET = 
  CPPUNIT_LDFLAGS_SET = 
endif
AM_CPPFLAGS += $(CPPUNIT_CPPFLAGS_SET)
AM_LDFLAGS += $(CPPUNIT_LDFLAGS_SET)

no_cpp_unit:
	@echo "echo '###################################################'"  > $@
	@echo "echo '#### Mesquite configured w/out CppUnit support ####'" >> $@
	@echo "echo '####          Cannot build unit tests          ####'" >> $@
	@echo "echo '###################################################'" >> $@
	@echo "exit 0" >> $@
	@chmod a+x $@

if ENABLE_IMESH
  IMESH_TEST_SRC = iMeshTest.cpp
  IMESH_LDADD = $(top_builddir)/itaps/libmsqimesh.la $(IMESH_LIBS)
  IMESH_DEP = $(top_builddir)/itaps/libmsqimesh.la
else
  IMESH_TEST_SRC =
  IMESH_LDADD = 
  IMESH_DEP =
endif
LDADD += $(IMESH_LDADD)
msq_test_DEPENDENCIES += $(IMESH_DEP)

if WITH_EXODUS
  EXODUS_TEST_SRC = ExodusTest.cpp
else
  EXODUS_TEST_SRC = 
endif
  
msq_test_SOURCES = msq_test_main.cpp \
	MesquiteTestRunner.hpp \
	ObjectiveFunctionTests.hpp \
	ObjectiveFunctionTests.cpp \
	PatchDataInstances.hpp \
	QualityMetricTester.hpp \
	QualityMetricTester.cpp \
	UnitUtil.hpp \
        TargetMetricTest.hpp \
	TMPQualityMetricTest.hpp \
	ArrayMeshTest.cpp \
	AspectRatioGammaTest.cpp \
	AveragingQMTest.cpp \
	AWMetricTest.cpp \
        AWQualityMetricTest.cpp \
	BCDTest.cpp \
	BoundedCylinderDomainTest.cpp \
	CachingTargetTest.cpp \
	CircleDomainTest.cpp \
	CompositeMetricTest.cpp \
	CompositeOFTest.cpp \
	ConditionNumberTest.cpp \
	ConicDomainTest.cpp \
	CylinderDomainTest.cpp \
        DomainClassifierTest.cpp \
	EdgeLengthMetricTest.cpp \
	EdgeLengthQualityMetricTest.cpp \
	ElemSampleQMTest.cpp \
	ExtraDataTest.cpp \
	FileTokenizerTest.cpp \
        GeomPrimTest.cpp \
	HexLagrangeShapeTest.cpp \
	HigherOrderTest.cpp \
	IdealElementTest.cpp \
	IdealTargetTest.cpp \
	IdealWeightMeanRatioTest.cpp \
	IdealWeightInverseMeanRatioTest.cpp \
	InstructionQueueTest.cpp \
	LinearMappingFunctionTest.cpp \
	LineDomainTest.cpp \
	LVQDTargetTest.cpp \
	MappingFunctionTest.cpp \
	Matrix3DTest.cpp \
	MeshInterfaceTest.cpp \
	MeshImplTest.cpp \
	MeshUtilTest.cpp \
	MiscTargetMetricTests.cpp \
	MsqFreeVertexIndexIteratorTest.cpp \
	MsqHessianTest.cpp \
	MsqMatrixTest.cpp \
	MsqMeshEntityTest.cpp \
	MsqVertexTest.cpp \
        NodeSetTest.cpp \
	NumericalOFTest.cpp \
	ObjectiveFunctionTest.cpp \
	PatchDataTest.cpp \
	PatchDataTestNormals.cpp \
	PatchSetTest.cpp \
	PlanarGeometryTest.cpp \
	PMeanPMetricTest.cpp \
	PMeanPTemplateTest.cpp \
	QuadLagrangeShapeTest.cpp \
	QualityAssessorTest.cpp \
	QualityMetricTest.cpp \
	ShapeOrientTargetMetricTests.cpp \
	ShapeSizeOrientTargetMetricTests.cpp \
	ShapeSizeTargetMetricTests.cpp \
	ShapeTargetMetricTests.cpp \
	SizeTargetMetricTests.cpp \
	SlaveBoundaryVerticesTest.cpp \
	SphericalDomainTest.cpp \
	SphericalGeometryTest.cpp \
	StdDevTemplateTest.cpp \
	SymMatrix3DTest.cpp \
	TagVertexMeshTest.cpp \
	TargetCalculatorTest.cpp \
	TargetReadWriteTest.cpp \
	TerminationCriterionTest.cpp \
	TetLagrangeShapeTest.cpp \
	TriLagrangeShapeTest.cpp \
	TopologyInfoTest.cpp \
	TMetricTest.cpp \
        TMPDerivsTest.cpp \
        TQualityMetricTest.cpp \
	UntangleBetaTest.cpp \
	UntangleTargetMetricTests.cpp \
	Vector3DTest.cpp \
	VertexConditionNumberTest.cpp \
	VertexCullingRegressionTest.cpp \
	VtkTest.cpp \
	XYRectangleTest.cpp \
        $(IMESH_TEST_SRC) \
        $(EXODUS_TEST_SRC)

CLEANFILES = no_cpp_unit
MOSTLYCLEANFILES = hexagon.vtk \
		   LPtoP-hex-bcd.vtk \
		   LPtoP-hex-gbl.vtk \
		   LPtoP-tet-bcd.vtk \
		   LPtoP-tet-gbl.vtk \
		   PMeanP-hex-bcd.vtk \
		   PMeanP-hex-gbl.vtk \
		   PMeanP-tet-bcd.vtk \
		   PMeanP-tet-gbl.vtk \
		   domain_ideal.vtk \
		   domain_ref.vtk \
		   none_ideal.vtk \
		   none_ref.vtk \
		   target_ideal.vtk \
		   target_ref.vtk

EXTRA_DIST = CMakeLists.txt

		   
