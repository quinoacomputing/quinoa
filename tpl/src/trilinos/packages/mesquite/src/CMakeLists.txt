IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")


  #
  # A) Package-specific configuration options
  #

ELSE()
  MACRO(APPEND_SET VARNAME)
    SET(${VARNAME} ${${VARNAME}} ${ARGN})
  ENDMACRO()
  
  FUNCTION(TRIBITS_ADD_LIBRARY name)
    LIST(REMOVE_ITEM ARGV "${name}" "HEADERS" "SOURCES")
    IF(${PACKAGE_NAME}_ENABLE_SHARED)
      ADD_LIBRARY( ${name} SHARED ${ARGV} )
    ELSE()
      ADD_LIBRARY( ${name} STATIC ${ARGV} )
    ENDIF()
  ENDFUNCTION()
ENDIF()



autoconf_header( ${CMAKE_CURRENT_SOURCE_DIR}/include/mesquite_version.h.in
                 ${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_version.h )
autoconf_header( ${CMAKE_CURRENT_SOURCE_DIR}/include/mesquite_config.h.in
                 ${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_config.h )

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_version.h
  include/Mesquite.hpp
  include/Mesquite_MeshInterface.hpp
  include/Mesquite_ParallelMeshInterface.hpp
  include/Mesquite_ParallelHelperInterface.hpp
  include/Mesquite_CurveDomain.hpp
  )

#
# src/Mesh
#

IF (TPL_ENABLE_MPI)
  APPEND_SET(HEADERS
    Mesh/Mesquite_ParallelHelper.hpp
    )

  APPEND_SET(SOURCES
    Mesh/Mesquite_ParallelHelper.cpp
    )
ELSE()
  APPEND_SET(HEADERS
    Mesh/Mesquite_ParallelHelperNoMPIStub.hpp
    )

  APPEND_SET(SOURCES
    Mesh/Mesquite_ParallelHelperNoMPIStub.cpp
    )
ENDIF()

APPEND_SET(HEADERS
  Mesh/Mesquite_ArrayMesh.hpp
  Mesh/Mesquite_ElementPatches.hpp
  Mesh/Mesquite_ExtraData.hpp
  Mesh/Mesquite_ExtraDataUser.hpp
  Mesh/Mesquite_GlobalPatch.hpp
  Mesh/Mesquite_IdealElements.hpp
  Mesh/Mesquite_MeshDecorator.hpp
  Mesh/Mesquite_MeshImpl.hpp
  Mesh/Mesquite_MeshImplData.hpp
  Mesh/Mesquite_MeshImplTags.hpp
  Mesh/Mesquite_MeshUtil.hpp
  Mesh/Mesquite_MeshWriter.hpp
  Mesh/Mesquite_MsqFreeVertexIndexIterator.hpp
  Mesh/Mesquite_MsqMeshEntity.hpp
  Mesh/Mesquite_MsqVertex.hpp
  Mesh/Mesquite_ParallelMeshImpl.hpp
  Mesh/Mesquite_PatchData.hpp
  Mesh/Mesquite_PatchIterator.hpp
  Mesh/Mesquite_PatchSet.hpp
  Mesh/Mesquite_TagVertexMesh.hpp
  Mesh/Mesquite_TopologyInfo.hpp
  Mesh/Mesquite_VertexPatches.hpp
  )

APPEND_SET(SOURCES
  Mesh/Mesquite_ArrayMesh.cpp
  Mesh/Mesquite_ElementPatches.cpp
  Mesh/Mesquite_ExtraData.cpp
  Mesh/Mesquite_GlobalPatch.cpp
  Mesh/Mesquite_IdealElements.cpp
  Mesh/Mesquite_MeshDecorator.cpp
  Mesh/Mesquite_MeshImpl.cpp
  Mesh/Mesquite_MeshImplData.cpp
  Mesh/Mesquite_MeshImplTags.cpp
  Mesh/Mesquite_MeshUtil.cpp
  Mesh/Mesquite_MeshWriter.cpp
  Mesh/Mesquite_MsqMeshEntity.cpp
  Mesh/Mesquite_MsqVertex.cpp
  Mesh/Mesquite_ParallelMeshImpl.cpp
  Mesh/Mesquite_PatchData.cpp
  Mesh/Mesquite_PatchIterator.cpp
  Mesh/Mesquite_PatchSet.cpp
  Mesh/Mesquite_TagVertexMesh.cpp
  Mesh/Mesquite_TopologyInfo.cpp
  Mesh/Mesquite_VertexPatches.cpp
  )

#
# src/Control
#

APPEND_SET(HEADERS
  Control/Mesquite_Instruction.hpp
  Control/Mesquite_InstructionQueue.hpp
  Control/Mesquite_IQInterface.hpp
  Control/Mesquite_TerminationCriterion.hpp
  Control/Mesquite_Settings.hpp
  Control/Mesquite_SlaveBoundaryVertices.hpp
  Control/Mesquite_VertexSlaver.hpp 
  )

APPEND_SET(SOURCES
  Control/Mesquite_Instruction.cpp
  Control/Mesquite_InstructionQueue.cpp
  Control/Mesquite_IQInterface.cpp
  Control/Mesquite_TerminationCriterion.cpp
  Control/Mesquite_Settings.cpp
  Control/Mesquite_SlaveBoundaryVertices.cpp 
  )


#
# src/Wrappers
#

APPEND_SET(HEADERS
  Wrappers/Mesquite_DeformingDomainWrapper.hpp
  Wrappers/Mesquite_LaplaceWrapper.hpp
  Wrappers/Mesquite_ShapeImprovementWrapper.hpp
  Wrappers/Mesquite_ShapeImprover.hpp
  Wrappers/Mesquite_SizeAdaptShapeWrapper.hpp
  Wrappers/Mesquite_PaverMinEdgeLengthWrapper.hpp 
  Wrappers/Mesquite_UntangleWrapper.hpp 
  Wrappers/Mesquite_ViscousCFDTetShapeWrapper.hpp 
  Wrappers/Mesquite_Wrapper.hpp 
  )

APPEND_SET(SOURCES
  Wrappers/Mesquite_DeformingDomainWrapper.cpp
  Wrappers/Mesquite_LaplaceWrapper.cpp
  Wrappers/Mesquite_ShapeImprovementWrapper.cpp
  Wrappers/Mesquite_ShapeImprover.cpp
  Wrappers/Mesquite_SizeAdaptShapeWrapper.cpp
  Wrappers/Mesquite_PaverMinEdgeLengthWrapper.cpp
  Wrappers/Mesquite_UntangleWrapper.cpp 
  Wrappers/Mesquite_ViscousCFDTetShapeWrapper.cpp
  Wrappers/Mesquite_Wrapper.cpp 
  )

#
# src/MappingFunction
#

APPEND_SET(HEADERS
  MappingFunction/Mesquite_MappingFunction.hpp
  MappingFunction/Mesquite_NodeSet.hpp
  MappingFunction/Mesquite_Sample.hpp
  )

APPEND_SET(SOURCES
  MappingFunction/Mesquite_MappingFunction.cpp
  MappingFunction/Mesquite_NodeSet.cpp
  MappingFunction/Mesquite_Sample.cpp
  )

#
# src/MappingFunction/Lagrange
#

APPEND_SET(HEADERS
  MappingFunction/Lagrange/Mesquite_HexLagrangeShape.hpp
  MappingFunction/Lagrange/Mesquite_QuadLagrangeShape.hpp
  MappingFunction/Lagrange/Mesquite_TetLagrangeShape.hpp
  MappingFunction/Lagrange/Mesquite_TriLagrangeShape.hpp 
  )

APPEND_SET(SOURCES
  MappingFunction/Lagrange/Mesquite_HexLagrangeShape.cpp
  MappingFunction/Lagrange/Mesquite_QuadLagrangeShape.cpp
  MappingFunction/Lagrange/Mesquite_TetLagrangeShape.cpp
  MappingFunction/Lagrange/Mesquite_TriLagrangeShape.cpp 
  )

#
# src/MappingFunction/Linear
#

APPEND_SET(HEADERS
  MappingFunction/Linear/Mesquite_LinearHexahedron.hpp
  MappingFunction/Linear/Mesquite_LinearPrism.hpp
  MappingFunction/Linear/Mesquite_LinearPyramid.hpp
  MappingFunction/Linear/Mesquite_LinearQuadrilateral.hpp
  MappingFunction/Linear/Mesquite_LinearTetrahedron.hpp
  MappingFunction/Linear/Mesquite_LinearTriangle.hpp 
  )

APPEND_SET(SOURCES
  MappingFunction/Linear/Mesquite_LinearHexahedron.cpp
  MappingFunction/Linear/Mesquite_LinearPrism.cpp
  MappingFunction/Linear/Mesquite_LinearPyramid.cpp
  MappingFunction/Linear/Mesquite_LinearQuadrilateral.cpp
  MappingFunction/Linear/Mesquite_LinearTetrahedron.cpp
  MappingFunction/Linear/Mesquite_LinearTriangle.cpp 
  )

#
# src/Misc
#

APPEND_SET(HEADERS
  Misc/Mesquite_BoundedCylinderDomain.hpp
  Misc/Mesquite_Bits.hpp
  Misc/Mesquite_ConicDomain.hpp
  Misc/Mesquite_CylinderDomain.hpp
  Misc/Mesquite_DomainClassifier.hpp
  Misc/Mesquite_EdgeIterator.hpp
  Misc/Mesquite_Exponent.hpp
  Misc/Mesquite_FileTokenizer.hpp
  Misc/Mesquite_Matrix3D.hpp
  Misc/Mesquite_MeshDomain1D.hpp
  Misc/Mesquite_MeshTransform.hpp
  Misc/Mesquite_MsqDebug.hpp
  Misc/Mesquite_MsqError.hpp
  Misc/Mesquite_MsqFPE.hpp
  Misc/Mesquite_MsqGeomPrim.hpp
  Misc/Mesquite_MsqHessian.hpp
  Misc/Mesquite_MsqInterrupt.hpp
  Misc/Mesquite_MsqMatrix.hpp
  Misc/Mesquite_MsqTimer.hpp
  Misc/Mesquite_PlanarDomain.hpp
  Misc/Mesquite_SimpleStats.hpp
  Misc/Mesquite_SphericalDomain.hpp
  Misc/Mesquite_SymMatrix3D.hpp
  Misc/Mesquite_Vector3D.hpp
  Misc/Mesquite_VtkTypeInfo.hpp
  Misc/Mesquite_XYPlanarDomain.hpp
  Misc/Mesquite_XYRectangle.hpp
  )

APPEND_SET(SOURCES
  Misc/Mesquite_BoundedCylinderDomain.cpp
  Misc/Mesquite_ConicDomain.cpp
  Misc/Mesquite_CylinderDomain.cpp
  Misc/Mesquite_DomainClassifier.cpp
  Misc/Mesquite_DomainUtil.cpp
  Misc/Mesquite_DomainUtil.hpp
  Misc/Mesquite_EdgeIterator.cpp
  Misc/Mesquite_Exponent.cpp
  Misc/Mesquite_FileTokenizer.cpp
  Misc/Mesquite_MeshDomain1D.cpp
  Misc/Mesquite_MeshTransform.cpp
  Misc/Mesquite_MesquiteVersion.cpp
  Misc/Mesquite_MsqDebug.cpp
  Misc/Mesquite_MsqError.cpp
  Misc/Mesquite_MsqFPE.cpp
  Misc/Mesquite_MsqGeomPrim.cpp
  Misc/Mesquite_MsqHessian.cpp
  Misc/Mesquite_MsqInterrupt.cpp
  Misc/Mesquite_MsqTimer.cpp
  Misc/Mesquite_PlanarDomain.cpp
  Misc/Mesquite_SimpleStats.cpp
  Misc/Mesquite_SphericalDomain.cpp
  Misc/Mesquite_Vector3D.cpp
  Misc/Mesquite_VtkTypeInfo.cpp
  Misc/Mesquite_XYPlanarDomain.cpp
  Misc/Mesquite_XYRectangle.cpp
  )
IF(IMESH_DEFS OR IGEOM_DEFS OR IREL_DEFS)
  APPEND_SET(SOURCES Misc/Mesquite_MsqIBase.cpp Misc/Mesquite_MsqIBase.hpp)
ENDIF()


#
# src/ObjectiveFuntion
#

APPEND_SET(HEADERS
  ObjectiveFunction/Mesquite_CompositeOFAdd.hpp
  ObjectiveFunction/Mesquite_CompositeOFMultiply.hpp
  ObjectiveFunction/Mesquite_CompositeOFScalarAdd.hpp
  ObjectiveFunction/Mesquite_CompositeOFScalarMultiply.hpp
  ObjectiveFunction/Mesquite_LInfTemplate.hpp
  ObjectiveFunction/Mesquite_LPtoPTemplate.hpp
  ObjectiveFunction/Mesquite_MaxTemplate.hpp
  ObjectiveFunction/Mesquite_OFEvaluator.hpp
  ObjectiveFunction/Mesquite_ObjectiveFunction.hpp
  ObjectiveFunction/Mesquite_ObjectiveFunctionTemplate.hpp
  ObjectiveFunction/Mesquite_PatchPowerMeanP.hpp
  ObjectiveFunction/Mesquite_PMeanPTemplate.hpp
  ObjectiveFunction/Mesquite_StdDevTemplate.hpp
  ObjectiveFunction/Mesquite_VarianceTemplate.hpp
  )

APPEND_SET(SOURCES
  ObjectiveFunction/Mesquite_CompositeOFAdd.cpp
  ObjectiveFunction/Mesquite_CompositeOFMultiply.cpp
  ObjectiveFunction/Mesquite_CompositeOFScalarAdd.cpp
  ObjectiveFunction/Mesquite_CompositeOFScalarMultiply.cpp
  ObjectiveFunction/Mesquite_LInfTemplate.cpp
  ObjectiveFunction/Mesquite_LPtoPTemplate.cpp
  ObjectiveFunction/Mesquite_MaxTemplate.cpp
  ObjectiveFunction/Mesquite_OFEvaluator.cpp
  ObjectiveFunction/Mesquite_ObjectiveFunction.cpp
  ObjectiveFunction/Mesquite_ObjectiveFunctionTemplate.cpp
  ObjectiveFunction/Mesquite_PatchPowerMeanP.cpp
  ObjectiveFunction/Mesquite_PMeanPTemplate.cpp
  ObjectiveFunction/Mesquite_StdDevTemplate.cpp
  ObjectiveFunction/Mesquite_VarianceTemplate.cpp
  )

#
# src/QualityAssessor
#

APPEND_SET(HEADERS
  QualityAssessor/Mesquite_QualityAssessor.hpp
  )

APPEND_SET(SOURCES
  QualityAssessor/Mesquite_QualityAssessor.cpp
  )

#
# src/QualityImprover
#

APPEND_SET(HEADERS
  QualityImprover/Mesquite_NullImprover.hpp
  QualityImprover/Mesquite_PatchSetUser.hpp
  QualityImprover/Mesquite_QualityImprover.hpp
  QualityImprover/Mesquite_VertexMover.hpp
  )

APPEND_SET(SOURCES
  QualityImprover/Mesquite_PatchSetUser.cpp
  QualityImprover/Mesquite_QualityImprover.cpp
  QualityImprover/Mesquite_VertexMover.cpp
  )

#
# src/QualityImprover/OptSolvers
#

APPEND_SET(HEADERS
  QualityImprover/OptSolvers/Mesquite_ConjugateGradient.hpp
  QualityImprover/OptSolvers/Mesquite_FeasibleNewton.hpp
  QualityImprover/OptSolvers/Mesquite_NonSmoothDescent.hpp
  QualityImprover/OptSolvers/Mesquite_QuasiNewton.hpp
  QualityImprover/OptSolvers/Mesquite_SteepestDescent.hpp
  QualityImprover/OptSolvers/Mesquite_NonGradient.hpp
  QualityImprover/OptSolvers/Mesquite_TrustRegion.hpp
  )

APPEND_SET(SOURCES
  QualityImprover/OptSolvers/Mesquite_ConjugateGradient.cpp
  QualityImprover/OptSolvers/Mesquite_FeasibleNewton.cpp
  QualityImprover/OptSolvers/Mesquite_NonSmoothDescent.cpp
  QualityImprover/OptSolvers/Mesquite_QuasiNewton.cpp
  QualityImprover/OptSolvers/Mesquite_SteepestDescent.cpp
  QualityImprover/OptSolvers/Mesquite_NonGradient.cpp
  QualityImprover/OptSolvers/Mesquite_TrustRegion.cpp
  )

#
# src/QualityImprover/Relaxation
#

APPEND_SET(HEADERS
  QualityImprover/Relaxation/Mesquite_Randomize.hpp
  QualityImprover/Relaxation/Mesquite_RelaxationSmoother.hpp
  QualityImprover/Relaxation/Mesquite_LaplacianSmoother.hpp
  QualityImprover/Relaxation/Mesquite_SmartLaplacianSmoother.hpp
  )

APPEND_SET(SOURCES
  QualityImprover/Relaxation/Mesquite_Randomize.cpp
  QualityImprover/Relaxation/Mesquite_RelaxationSmoother.cpp
  QualityImprover/Relaxation/Mesquite_LaplacianSmoother.cpp
  QualityImprover/Relaxation/Mesquite_SmartLaplacianSmoother.cpp
  )

#
# src/QualityMetric
#

APPEND_SET(HEADERS
  QualityMetric/Mesquite_AddQualityMetric.hpp
  QualityMetric/Mesquite_AveragingQM.hpp
  QualityMetric/Mesquite_EdgeQM.hpp
  QualityMetric/Mesquite_ElemSampleQM.hpp
  QualityMetric/Mesquite_ElementMaxQM.hpp
  QualityMetric/Mesquite_ElementAvgQM.hpp
  QualityMetric/Mesquite_ElementPMeanP.hpp
  QualityMetric/Mesquite_ElementQM.hpp
  QualityMetric/Mesquite_MultiplyQualityMetric.hpp
  QualityMetric/Debug/Mesquite_NumericalQM.hpp
  QualityMetric/Debug/Mesquite_CompareQM.hpp
  QualityMetric/Mesquite_PMeanPMetric.hpp
  QualityMetric/Mesquite_PowerQualityMetric.hpp
  QualityMetric/Mesquite_QualityMetric.hpp
  QualityMetric/Mesquite_ScalarAddQualityMetric.hpp
  QualityMetric/Mesquite_ScalarMultiplyQualityMetric.hpp
  QualityMetric/Mesquite_VertexMaxQM.hpp
  QualityMetric/Mesquite_VertexPMeanP.hpp
  QualityMetric/Mesquite_VertexQM.hpp 
  )

APPEND_SET(SOURCES
  QualityMetric/Mesquite_AddQualityMetric.cpp
  QualityMetric/Mesquite_AveragingQM.cpp
  QualityMetric/Mesquite_EdgeQM.cpp
  QualityMetric/Mesquite_ElemSampleQM.cpp
  QualityMetric/Mesquite_ElementMaxQM.cpp
  QualityMetric/Mesquite_ElementAvgQM.cpp
  QualityMetric/Mesquite_ElementPMeanP.cpp
  QualityMetric/Mesquite_ElementQM.cpp
  QualityMetric/Mesquite_MultiplyQualityMetric.cpp
  QualityMetric/Debug/Mesquite_NumericalQM.cpp
  QualityMetric/Debug/Mesquite_CompareQM.cpp
  QualityMetric/Mesquite_PMeanPMetric.cpp
  QualityMetric/Mesquite_PowerQualityMetric.cpp
  QualityMetric/Mesquite_QualityMetric.cpp
  QualityMetric/Mesquite_ScalarAddQualityMetric.cpp
  QualityMetric/Mesquite_ScalarMultiplyQualityMetric.cpp
  QualityMetric/Mesquite_VertexMaxQM.cpp
  QualityMetric/Mesquite_VertexPMeanP.cpp
  QualityMetric/Mesquite_VertexQM.cpp 
  )

#
# src/QualityMetric/Shape
#

APPEND_SET(HEADERS
  QualityMetric/Shape/Mesquite_AspectRatioGammaQualityMetric.hpp
  QualityMetric/Shape/Mesquite_ConditionNumberFunctions.hpp
  QualityMetric/Shape/Mesquite_ConditionNumberQualityMetric.hpp
  QualityMetric/Shape/Mesquite_IdealWeightInverseMeanRatio.hpp
  QualityMetric/Shape/Mesquite_IdealWeightMeanRatio.hpp
  QualityMetric/Shape/Mesquite_MeanRatioFunctions.hpp
  QualityMetric/Shape/Mesquite_VertexConditionNumberQualityMetric.hpp
  )

APPEND_SET(SOURCES
  QualityMetric/Shape/Mesquite_AspectRatioGammaQualityMetric.cpp
  QualityMetric/Shape/Mesquite_ConditionNumberQualityMetric.cpp
  QualityMetric/Shape/Mesquite_IdealWeightInverseMeanRatio.cpp
  QualityMetric/Shape/Mesquite_IdealWeightMeanRatio.cpp
  QualityMetric/Shape/Mesquite_VertexConditionNumberQualityMetric.cpp
  )

#
# src/QualityMetric/Smoothness
#

APPEND_SET(HEADERS
  QualityMetric/Smoothness/Mesquite_EdgeLengthQualityMetric.hpp
  QualityMetric/Smoothness/Mesquite_EdgeLengthRangeQualityMetric.hpp
  )

APPEND_SET(SOURCES
  QualityMetric/Smoothness/Mesquite_EdgeLengthQualityMetric.cpp
  QualityMetric/Smoothness/Mesquite_EdgeLengthRangeQualityMetric.cpp
  )

#
# src/QualityMetric/TMP
#

APPEND_SET(HEADERS
  QualityMetric/TMP/Mesquite_AffineMapMetric.hpp
  QualityMetric/TMP/Mesquite_TMPQualityMetric.hpp
  QualityMetric/TMP/Mesquite_AWQualityMetric.hpp 
  QualityMetric/TMP/Mesquite_TQualityMetric.hpp 
  )

APPEND_SET(SOURCES
  QualityMetric/TMP/Mesquite_AffineMapMetric.cpp
  QualityMetric/TMP/Mesquite_TargetMetricUtil.cpp
  QualityMetric/TMP/Mesquite_TargetMetricUtil.hpp
  QualityMetric/TMP/Mesquite_TMPQualityMetric.cpp 
  QualityMetric/TMP/Mesquite_AWQualityMetric.cpp 
  QualityMetric/TMP/Mesquite_TQualityMetric.cpp 
  )

#
# src/QualityMetric/Untangle
#

APPEND_SET(HEADERS
  QualityMetric/Untangle/Mesquite_UntangleBetaQualityMetric.hpp
  )

APPEND_SET(SOURCES
  QualityMetric/Untangle/Mesquite_UntangleBetaQualityMetric.cpp
  )

#
# src/QualityMetric/Volume
#

APPEND_SET(HEADERS
  QualityMetric/Volume/Mesquite_EdgeLengthMetric.hpp
  QualityMetric/Volume/Mesquite_LocalSizeQualityMetric.hpp
  QualityMetric/Volume/Mesquite_SizeMetric.hpp
  )

APPEND_SET(SOURCES
  QualityMetric/Volume/Mesquite_EdgeLengthMetric.cpp
  QualityMetric/Volume/Mesquite_LocalSizeQualityMetric.cpp
  QualityMetric/Volume/Mesquite_SizeMetric.cpp
  )


#
# src/TargetMetric
#

APPEND_SET(HEADERS
  TargetMetric/Mesquite_AWMetric.hpp
  TargetMetric/Mesquite_TMetric.hpp
  TargetMetric/Mesquite_TMetricBarrier.hpp
  TargetMetric/Mesquite_TMetricNonBarrier.hpp
  TargetMetric/Mesquite_AWMetricBarrier.hpp
  TargetMetric/Mesquite_AWMetricNonBarrier.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/Mesquite_TMPCommon.hpp
  TargetMetric/Mesquite_TMPDerivs.hpp
  TargetMetric/Mesquite_AWMetric.cpp
  TargetMetric/Mesquite_TMetric.cpp
  TargetMetric/Mesquite_TMetricBarrier.cpp
  TargetMetric/Mesquite_TMetricNonBarrier.cpp
  TargetMetric/Mesquite_AWMetricBarrier.cpp
  TargetMetric/Mesquite_AWMetricNonBarrier.cpp
  )


#
# src/TargetMetric/Misc
#

APPEND_SET(HEADERS
  TargetMetric/Misc/Mesquite_InvTransBarrier.hpp
  TargetMetric/Misc/Mesquite_TMixed.hpp
  TargetMetric/Misc/Mesquite_TOffset.hpp
  TargetMetric/Misc/Mesquite_TPower2.hpp
  TargetMetric/Misc/Mesquite_TScale.hpp
  TargetMetric/Misc/Mesquite_TSquared.hpp 
  TargetMetric/Misc/Mesquite_TSum.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/Misc/Mesquite_InvTransBarrier.cpp
  TargetMetric/Misc/Mesquite_TMixed.cpp
  TargetMetric/Misc/Mesquite_TOffset.cpp
  TargetMetric/Misc/Mesquite_TPower2.cpp
  TargetMetric/Misc/Mesquite_TScale.cpp
  TargetMetric/Misc/Mesquite_TSquared.cpp 
  TargetMetric/Misc/Mesquite_TSum.cpp
  )

#
# src/TargetMetric/Shape
#

APPEND_SET(HEADERS
  TargetMetric/Shape/Mesquite_AWShape2DB1.hpp
  TargetMetric/Shape/Mesquite_AWShape2DNB1.hpp
  TargetMetric/Shape/Mesquite_AWShape2DNB2.hpp
  TargetMetric/Shape/Mesquite_TInverseMeanRatio.hpp
  TargetMetric/Shape/Mesquite_TShape2DNB2.hpp
  TargetMetric/Shape/Mesquite_TShape3DB2.hpp
  TargetMetric/Shape/Mesquite_TShapeB1.hpp
  TargetMetric/Shape/Mesquite_TShapeNB1.hpp 
  )

APPEND_SET(SOURCES
  TargetMetric/Shape/Mesquite_AWShape2DB1.cpp
  TargetMetric/Shape/Mesquite_AWShape2DNB1.cpp
  TargetMetric/Shape/Mesquite_AWShape2DNB2.cpp
  TargetMetric/Shape/Mesquite_TInverseMeanRatio.cpp
  TargetMetric/Shape/Mesquite_TShape2DNB2.cpp
  TargetMetric/Shape/Mesquite_TShape3DB2.cpp
  TargetMetric/Shape/Mesquite_TShapeB1.cpp
  TargetMetric/Shape/Mesquite_TShapeNB1.cpp 
  )

#
# src/TargetMetric/ShapeOrient
#

APPEND_SET(HEADERS
  TargetMetric/ShapeOrient/Mesquite_AWShapeOrientNB1.hpp 
  TargetMetric/ShapeOrient/Mesquite_TShapeOrientB1.hpp 
  TargetMetric/ShapeOrient/Mesquite_TShapeOrientB2.hpp 
  TargetMetric/ShapeOrient/Mesquite_TShapeOrientNB1.hpp
  TargetMetric/ShapeOrient/Mesquite_TShapeOrientNB2.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/ShapeOrient/Mesquite_AWShapeOrientNB1.cpp 
  TargetMetric/ShapeOrient/Mesquite_TShapeOrientB1.cpp 
  TargetMetric/ShapeOrient/Mesquite_TShapeOrientB2.cpp 
  TargetMetric/ShapeOrient/Mesquite_TShapeOrientNB1.cpp
  TargetMetric/ShapeOrient/Mesquite_TShapeOrientNB2.cpp
  )

#
# src/TargetMetric/ShapeSize
#

APPEND_SET(HEADERS
  TargetMetric/ShapeSize/Mesquite_AWShapeSizeB1.hpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize2DB2.hpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize2DNB1.hpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize2DNB2.hpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize3DB2.hpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize3DB4.hpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize3DNB1.hpp
  TargetMetric/ShapeSize/Mesquite_TShapeSizeB1.hpp
  TargetMetric/ShapeSize/Mesquite_TShapeSizeB3.hpp
  TargetMetric/ShapeSize/Mesquite_TShapeSizeNB3.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/ShapeSize/Mesquite_AWShapeSizeB1.cpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize2DB2.cpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize2DNB1.cpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize2DNB2.cpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize3DB2.cpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize3DB4.cpp
  TargetMetric/ShapeSize/Mesquite_TShapeSize3DNB1.cpp
  TargetMetric/ShapeSize/Mesquite_TShapeSizeB1.cpp
  TargetMetric/ShapeSize/Mesquite_TShapeSizeB3.cpp
  TargetMetric/ShapeSize/Mesquite_TShapeSizeNB3.cpp
  )

#
# src/TargetMetric/ShapeSizeOrient
#

APPEND_SET(HEADERS
  TargetMetric/ShapeSizeOrient/Mesquite_AWShapeSizeOrientNB1.hpp
  TargetMetric/ShapeSizeOrient/Mesquite_TShapeSizeOrientB1.hpp
  TargetMetric/ShapeSizeOrient/Mesquite_TShapeSizeOrientB2.hpp
  TargetMetric/ShapeSizeOrient/Mesquite_TShapeSizeOrientNB1.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/ShapeSizeOrient/Mesquite_AWShapeSizeOrientNB1.cpp
  TargetMetric/ShapeSizeOrient/Mesquite_TShapeSizeOrientB1.cpp
  TargetMetric/ShapeSizeOrient/Mesquite_TShapeSizeOrientB2.cpp
  TargetMetric/ShapeSizeOrient/Mesquite_TShapeSizeOrientNB1.cpp
  )

#
# src/TargetMetric/Size
#

APPEND_SET(HEADERS
  TargetMetric/Size/Mesquite_AWSizeB1.hpp
  TargetMetric/Size/Mesquite_AWSizeNB1.hpp
  TargetMetric/Size/Mesquite_TSizeB1.hpp
  TargetMetric/Size/Mesquite_TSizeNB1.hpp
  TargetMetric/Size/Mesquite_TTau.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/Size/Mesquite_AWSizeB1.cpp
  TargetMetric/Size/Mesquite_AWSizeNB1.cpp
  TargetMetric/Size/Mesquite_TSizeB1.cpp
  TargetMetric/Size/Mesquite_TSizeNB1.cpp
  TargetMetric/Size/Mesquite_TTau.cpp
  )

#
# src/TargetMetric/Untangle
#

APPEND_SET(HEADERS
  TargetMetric/Untangle/Mesquite_AWUntangleBeta.hpp
  TargetMetric/Untangle/Mesquite_TUntangle1.hpp
  TargetMetric/Untangle/Mesquite_TUntangleBeta.hpp
  TargetMetric/Untangle/Mesquite_TUntangleMu.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/Untangle/Mesquite_AWUntangleBeta.cpp
  TargetMetric/Untangle/Mesquite_TUntangle1.cpp
  TargetMetric/Untangle/Mesquite_TUntangleBeta.cpp
  TargetMetric/Untangle/Mesquite_TUntangleMu.cpp
  )

#
# src/TargetCalculator
#

APPEND_SET(HEADERS
  TargetCalculator/Mesquite_CachingTargetCalculator.hpp
  TargetCalculator/Mesquite_IdealShapeTarget.hpp
  TargetCalculator/Mesquite_InverseMetricWeight.hpp
  TargetCalculator/Mesquite_JacobianCalculator.hpp
  TargetCalculator/Mesquite_LambdaTarget.hpp
  TargetCalculator/Mesquite_LambdaConstant.hpp
  TargetCalculator/Mesquite_LVQDTargetCalculator.hpp
  TargetCalculator/Mesquite_MetricWeight.hpp
  TargetCalculator/Mesquite_RefMeshTargetCalculator.hpp
  TargetCalculator/Mesquite_ReferenceMesh.hpp
  TargetCalculator/Mesquite_RefSizeTargetCalculator.hpp
  TargetCalculator/Mesquite_RemainingWeight.hpp
  TargetCalculator/Mesquite_TargetCalculator.hpp
  TargetCalculator/Mesquite_TargetReader.hpp
  TargetCalculator/Mesquite_TargetWriter.hpp
  TargetCalculator/Mesquite_TetDihedralWeight.hpp
  TargetCalculator/Mesquite_WeightCalculator.hpp
  TargetCalculator/Mesquite_WeightReader.hpp 
  )

APPEND_SET(SOURCES
  TargetCalculator/Mesquite_CachingTargetCalculator.cpp
  TargetCalculator/Mesquite_IdealShapeTarget.cpp
  TargetCalculator/Mesquite_InverseMetricWeight.cpp
  TargetCalculator/Mesquite_JacobianCalculator.cpp
  TargetCalculator/Mesquite_LambdaTarget.cpp
  TargetCalculator/Mesquite_LambdaConstant.cpp
  TargetCalculator/Mesquite_LVQDTargetCalculator.cpp
  TargetCalculator/Mesquite_MetricWeight.cpp
  TargetCalculator/Mesquite_RefMeshTargetCalculator.cpp
  TargetCalculator/Mesquite_ReferenceMesh.cpp
  TargetCalculator/Mesquite_RefSizeTargetCalculator.cpp
  TargetCalculator/Mesquite_RemainingWeight.cpp
  TargetCalculator/Mesquite_TargetCalculator.cpp
  TargetCalculator/Mesquite_TargetReader.cpp
  TargetCalculator/Mesquite_TargetWriter.cpp
  TargetCalculator/Mesquite_TetDihedralWeight.cpp
  TargetCalculator/Mesquite_WeightCalculator.cpp
  TargetCalculator/Mesquite_WeightReader.cpp 
  )

#
# itaps
#


#the key here is that the generated files have to be marked as generated using the SET_SOURCE_FILES_PROPERTIES
#command. 
#
# Generate Fortran name mangling headers
#
# ADD_CUSTOM_COMMAND( 
#   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMesh_protos.h
#   COMMAND ${CMAKE_COMMAND}
#   ARGS    -D input_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/../itaps/iMesh.h 
#           -D output_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMesh_protos.h 
#           -D prefix:STRING=iMesh 
#           -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/ITAPSFortranMangling.cmake
#   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../itaps/iMesh.h
#   COMMENT "Running custom command"
#   )
# SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMesh_protos.h PROPERTIES GENERATED true)
# 
# ADD_CUSTOM_COMMAND( 
#   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMeshP_protos.h
#   COMMAND ${CMAKE_COMMAND}
#   ARGS    -D input_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/../itaps/iMeshP.h 
#           -D output_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMeshP_protos.h 
#           -D prefix:STRING=iMeshP 
#           -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/ITAPSFortranMangling.cmake
#   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../itaps/iMeshP.h
#   COMMENT "Running custom command"
#   )
# SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMeshP_protos.h PROPERTIES GENERATED true)
# 
# ADD_CUSTOM_COMMAND( 
#   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/../itaps/iGeom_protos.h
#   COMMAND ${CMAKE_COMMAND}
#   ARGS    -D input_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/../itaps/iGeom.h 
#           -D output_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/../itaps/iGeom_protos.h 
#           -D prefix:STRING=iGeom 
#           -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/ITAPSFortranMangling.cmake
#   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../itaps/iGeom.h
#   COMMENT "Running custom command"
#   )
# SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMeshP_protos.h PROPERTIES GENERATED true)
# 
# ADD_CUSTOM_COMMAND( 
#   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/../itaps/iRel_protos.h
#   COMMAND ${CMAKE_COMMAND}
#   ARGS    -D input_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/../itaps/iRel.h 
#           -D output_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/../itaps/iRel_protos.h 
#           -D prefix:STRING=iRel 
#           -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/ITAPSFortranMangling.cmake
#   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../itaps/iRel.h
#   COMMENT "Running custom command"
#   )
# SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMeshP_protos.h PROPERTIES GENERATED true)
# 
# 
#APPEND_SET(HEADERS
#  ../itaps/iBase.h
#  ../itaps/iBase_FCDefs.h
#  ../itaps/iMesh.h
#  ${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMesh_protos.h
#  ../itaps/iMeshP.h
#  ${CMAKE_CURRENT_BINARY_DIR}/../itaps/iMeshP_protos.h
#  ../itaps/iGeom.h
#  ${CMAKE_CURRENT_BINARY_DIR}/../itaps/iGeom_protos.h
#  ../itaps/iRel.h
#  ${CMAKE_CURRENT_BINARY_DIR}/../itaps/iRel_protos.h 
#  )


#
# Create the Mesquite_all_headers.hpp file.
#

SET( all_file ${CMAKE_CURRENT_BINARY_DIR}/include/Mesquite_all_headers.hpp )
# Cannot do this correctly because list if all headers is too long
# of a command-line argument for some platforms (e.g. Windows.)
# Instead, use less optimial solution (doesn't automatically update.)
#
#ADD_CUSTOM_COMMAND( 
#  OUTPUT ${all_file}
#  COMMAND ${CMAKE_COMMAND}
#  ARGS -D output_file:STRING=${all_file}
#       -D headers:STRING="${HEADERS}"
#       -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/GenerateAllHeaders.cmake
#  DEPENDS ${HEADERS}
#  COMMENT "Generating Mesquite_all_headers.hpp" )
#
SET_SOURCE_FILES_PROPERTIES(${all_file} PROPERTIES GENERATED true)
#
file( WRITE ${all_file} "#ifndef MESQUITE_ALL_HEADERS_HPP\n#define MESQUITE_ALL_HEADERS_HPP\n" )
foreach( header ${HEADERS} )
  string( REGEX MATCH "[^/]+\\.h(pp)?" file ${header} )
  file( APPEND ${all_file} "#include \"${file}\"\n")
endforeach()
file( APPEND ${all_file} "#endif\n" )
  


#last thing to do is add then newly created header
APPEND_SET(HEADERS ${all_file} )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  mesquite
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

MACRO(APPEND_DEFINE VARNAME OPT_FLAG)
  IF (${OPT_FLAG})
    SET(${VARNAME} "${${VARNAME}} -D${OPT_FLAG}")
  ENDIF()
ENDMACRO()

# s/MESQUITE_DLL_EXPORTS/MESQUITE_STATIC_LIB/ if building static library
IF(${PACKAGE_NAME}_ENABLE_SHARED)
  SET( MESQUITE_LIB_DEFINES "${MESQUITE_DEFINES} -DMESQUITE_DLL_EXPORTS" )
ELSE()
  SET( MESQUITE_LIB_DEFINES "${MESQUITE_DEFINES} -DMESQUITE_STATIC_LIB" )
ENDIF()
SET_SOURCE_FILES_PROPERTIES( ${SOURCES} COMPILE_FLAGS ${MESQUITE_LIB_DEFINES} )

