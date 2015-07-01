IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")

  INCLUDE(TribitsLibraryMacros)

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
  include/MeshInterface.hpp
  include/ParallelMeshInterface.hpp
  include/ParallelHelperInterface.hpp
  include/CurveDomain.hpp
  )

#
# src/Mesh
#

IF (TPL_ENABLE_MPI)
  APPEND_SET(HEADERS
    Mesh/ParallelHelper.hpp
    )

  APPEND_SET(SOURCES
    Mesh/ParallelHelper.cpp
    )
ELSE()
  APPEND_SET(HEADERS
    Mesh/ParallelHelperNoMPIStub.hpp
    )

  APPEND_SET(SOURCES
    Mesh/ParallelHelperNoMPIStub.cpp
    )
ENDIF()

APPEND_SET(HEADERS
  Mesh/ArrayMesh.hpp
  Mesh/ElementPatches.hpp
  Mesh/ExtraData.hpp
  Mesh/ExtraDataUser.hpp
  Mesh/GlobalPatch.hpp
  Mesh/IdealElements.hpp
  Mesh/MeshDecorator.hpp
  Mesh/MeshImpl.hpp
  Mesh/MeshImplData.hpp
  Mesh/MeshImplTags.hpp
  Mesh/MeshUtil.hpp
  Mesh/MeshWriter.hpp
  Mesh/MsqFreeVertexIndexIterator.hpp
  Mesh/MsqMeshEntity.hpp
  Mesh/MsqVertex.hpp
  Mesh/ParallelMeshImpl.hpp
  Mesh/PatchData.hpp
  Mesh/PatchIterator.hpp
  Mesh/PatchSet.hpp
  Mesh/TagVertexMesh.hpp
  Mesh/TopologyInfo.hpp
  Mesh/VertexPatches.hpp
  )

APPEND_SET(SOURCES
  Mesh/ArrayMesh.cpp
  Mesh/ElementPatches.cpp
  Mesh/ExtraData.cpp
  Mesh/GlobalPatch.cpp
  Mesh/IdealElements.cpp
  Mesh/MeshDecorator.cpp
  Mesh/MeshImpl.cpp
  Mesh/MeshImplData.cpp
  Mesh/MeshImplTags.cpp
  Mesh/MeshUtil.cpp
  Mesh/MeshWriter.cpp
  Mesh/MsqMeshEntity.cpp
  Mesh/MsqVertex.cpp
  Mesh/ParallelMeshImpl.cpp
  Mesh/PatchData.cpp
  Mesh/PatchIterator.cpp
  Mesh/PatchSet.cpp
  Mesh/TagVertexMesh.cpp
  Mesh/TopologyInfo.cpp
  Mesh/VertexPatches.cpp
  )

#
# src/Control
#

APPEND_SET(HEADERS
  Control/Instruction.hpp
  Control/InstructionQueue.hpp
  Control/IQInterface.hpp
  Control/TerminationCriterion.hpp
  Control/Settings.hpp
  Control/SlaveBoundaryVertices.hpp
  Control/VertexSlaver.hpp 
  )

APPEND_SET(SOURCES
  Control/Instruction.cpp
  Control/InstructionQueue.cpp
  Control/IQInterface.cpp
  Control/TerminationCriterion.cpp
  Control/Settings.cpp
  Control/SlaveBoundaryVertices.cpp 
  )


#
# src/Wrappers
#

APPEND_SET(HEADERS
  Wrappers/DeformingDomainWrapper.hpp
  Wrappers/LaplaceWrapper.hpp
  Wrappers/ShapeImprovementWrapper.hpp
  Wrappers/ShapeImprover.hpp
  Wrappers/SizeAdaptShapeWrapper.hpp
  Wrappers/PaverMinEdgeLengthWrapper.hpp 
  Wrappers/UntangleWrapper.hpp 
  Wrappers/ViscousCFDTetShapeWrapper.hpp 
  Wrappers/Wrapper.hpp 
  )

APPEND_SET(SOURCES
  Wrappers/DeformingDomainWrapper.cpp
  Wrappers/LaplaceWrapper.cpp
  Wrappers/ShapeImprovementWrapper.cpp
  Wrappers/ShapeImprover.cpp
  Wrappers/SizeAdaptShapeWrapper.cpp
  Wrappers/PaverMinEdgeLengthWrapper.cpp
  Wrappers/UntangleWrapper.cpp 
  Wrappers/ViscousCFDTetShapeWrapper.cpp
  Wrappers/Wrapper.cpp 
  )

#
# src/MappingFunction
#

APPEND_SET(HEADERS
  MappingFunction/MappingFunction.hpp
  MappingFunction/NodeSet.hpp
  MappingFunction/Sample.hpp
  )

APPEND_SET(SOURCES
  MappingFunction/MappingFunction.cpp
  MappingFunction/NodeSet.cpp
  MappingFunction/Sample.cpp
  )

#
# src/MappingFunction/Lagrange
#

APPEND_SET(HEADERS
  MappingFunction/Lagrange/HexLagrangeShape.hpp
  MappingFunction/Lagrange/QuadLagrangeShape.hpp
  MappingFunction/Lagrange/TetLagrangeShape.hpp
  MappingFunction/Lagrange/TriLagrangeShape.hpp 
  )

APPEND_SET(SOURCES
  MappingFunction/Lagrange/HexLagrangeShape.cpp
  MappingFunction/Lagrange/QuadLagrangeShape.cpp
  MappingFunction/Lagrange/TetLagrangeShape.cpp
  MappingFunction/Lagrange/TriLagrangeShape.cpp 
  )

#
# src/MappingFunction/Linear
#

APPEND_SET(HEADERS
  MappingFunction/Linear/LinearHexahedron.hpp
  MappingFunction/Linear/LinearPrism.hpp
  MappingFunction/Linear/LinearPyramid.hpp
  MappingFunction/Linear/LinearQuadrilateral.hpp
  MappingFunction/Linear/LinearTetrahedron.hpp
  MappingFunction/Linear/LinearTriangle.hpp 
  )

APPEND_SET(SOURCES
  MappingFunction/Linear/LinearHexahedron.cpp
  MappingFunction/Linear/LinearPrism.cpp
  MappingFunction/Linear/LinearPyramid.cpp
  MappingFunction/Linear/LinearQuadrilateral.cpp
  MappingFunction/Linear/LinearTetrahedron.cpp
  MappingFunction/Linear/LinearTriangle.cpp 
  )

#
# src/Misc
#

APPEND_SET(HEADERS
  Misc/BoundedCylinderDomain.hpp
  Misc/Bits.hpp
  Misc/ConicDomain.hpp
  Misc/CylinderDomain.hpp
  Misc/DomainClassifier.hpp
  Misc/EdgeIterator.hpp
  Misc/Exponent.hpp
  Misc/FileTokenizer.hpp
  Misc/Matrix3D.hpp
  Misc/MeshDomain1D.hpp
  Misc/MeshTransform.hpp
  Misc/MsqDebug.hpp
  Misc/MsqError.hpp
  Misc/MsqFPE.hpp
  Misc/MsqGeomPrim.hpp
  Misc/MsqHessian.hpp
  Misc/MsqInterrupt.hpp
  Misc/MsqMatrix.hpp
  Misc/MsqTimer.hpp
  Misc/PlanarDomain.hpp
  Misc/SimpleStats.hpp
  Misc/SphericalDomain.hpp
  Misc/SymMatrix3D.hpp
  Misc/Vector3D.hpp
  Misc/VtkTypeInfo.hpp
  Misc/XYPlanarDomain.hpp
  Misc/XYRectangle.hpp
  )

APPEND_SET(SOURCES
  Misc/BoundedCylinderDomain.cpp
  Misc/ConicDomain.cpp
  Misc/CylinderDomain.cpp
  Misc/DomainClassifier.cpp
  Misc/DomainUtil.cpp
  Misc/DomainUtil.hpp
  Misc/EdgeIterator.cpp
  Misc/Exponent.cpp
  Misc/FileTokenizer.cpp
  Misc/MeshDomain1D.cpp
  Misc/MeshTransform.cpp
  Misc/MesquiteVersion.cpp
  Misc/MsqDebug.cpp
  Misc/MsqError.cpp
  Misc/MsqFPE.cpp
  Misc/MsqGeomPrim.cpp
  Misc/MsqHessian.cpp
  Misc/MsqInterrupt.cpp
  Misc/MsqTimer.cpp
  Misc/PlanarDomain.cpp
  Misc/SimpleStats.cpp
  Misc/SphericalDomain.cpp
  Misc/Vector3D.cpp
  Misc/VtkTypeInfo.cpp
  Misc/XYPlanarDomain.cpp
  Misc/XYRectangle.cpp
  )
IF(IMESH_DEFS OR IGEOM_DEFS OR IREL_DEFS)
  APPEND_SET(SOURCES Misc/MsqIBase.cpp Misc/MsqIBase.hpp)
ENDIF()


#
# src/ObjectiveFuntion
#

APPEND_SET(HEADERS
  ObjectiveFunction/CompositeOFAdd.hpp
  ObjectiveFunction/CompositeOFMultiply.hpp
  ObjectiveFunction/CompositeOFScalarAdd.hpp
  ObjectiveFunction/CompositeOFScalarMultiply.hpp
  ObjectiveFunction/LInfTemplate.hpp
  ObjectiveFunction/LPtoPTemplate.hpp
  ObjectiveFunction/MaxTemplate.hpp
  ObjectiveFunction/OFEvaluator.hpp
  ObjectiveFunction/ObjectiveFunction.hpp
  ObjectiveFunction/ObjectiveFunctionTemplate.hpp
  ObjectiveFunction/PatchPowerMeanP.hpp
  ObjectiveFunction/PMeanPTemplate.hpp
  ObjectiveFunction/StdDevTemplate.hpp
  ObjectiveFunction/VarianceTemplate.hpp
  )

APPEND_SET(SOURCES
  ObjectiveFunction/CompositeOFAdd.cpp
  ObjectiveFunction/CompositeOFMultiply.cpp
  ObjectiveFunction/CompositeOFScalarAdd.cpp
  ObjectiveFunction/CompositeOFScalarMultiply.cpp
  ObjectiveFunction/LInfTemplate.cpp
  ObjectiveFunction/LPtoPTemplate.cpp
  ObjectiveFunction/MaxTemplate.cpp
  ObjectiveFunction/OFEvaluator.cpp
  ObjectiveFunction/ObjectiveFunction.cpp
  ObjectiveFunction/ObjectiveFunctionTemplate.cpp
  ObjectiveFunction/PatchPowerMeanP.cpp
  ObjectiveFunction/PMeanPTemplate.cpp
  ObjectiveFunction/StdDevTemplate.cpp
  ObjectiveFunction/VarianceTemplate.cpp
  )

#
# src/QualityAssessor
#

APPEND_SET(HEADERS
  QualityAssessor/QualityAssessor.hpp
  )

APPEND_SET(SOURCES
  QualityAssessor/QualityAssessor.cpp
  )

#
# src/QualityImprover
#

APPEND_SET(HEADERS
  QualityImprover/NullImprover.hpp
  QualityImprover/PatchSetUser.hpp
  QualityImprover/QualityImprover.hpp
  QualityImprover/VertexMover.hpp
  )

APPEND_SET(SOURCES
  QualityImprover/PatchSetUser.cpp
  QualityImprover/QualityImprover.cpp
  QualityImprover/VertexMover.cpp
  )

#
# src/QualityImprover/OptSolvers
#

APPEND_SET(HEADERS
  QualityImprover/OptSolvers/ConjugateGradient.hpp
  QualityImprover/OptSolvers/FeasibleNewton.hpp
  QualityImprover/OptSolvers/NonSmoothDescent.hpp
  QualityImprover/OptSolvers/QuasiNewton.hpp
  QualityImprover/OptSolvers/SteepestDescent.hpp
  QualityImprover/OptSolvers/NonGradient.hpp
  QualityImprover/OptSolvers/TrustRegion.hpp
  )

APPEND_SET(SOURCES
  QualityImprover/OptSolvers/ConjugateGradient.cpp
  QualityImprover/OptSolvers/FeasibleNewton.cpp
  QualityImprover/OptSolvers/NonSmoothDescent.cpp
  QualityImprover/OptSolvers/QuasiNewton.cpp
  QualityImprover/OptSolvers/SteepestDescent.cpp
  QualityImprover/OptSolvers/NonGradient.cpp
  QualityImprover/OptSolvers/TrustRegion.cpp
  )

#
# src/QualityImprover/Relaxation
#

APPEND_SET(HEADERS
  QualityImprover/Relaxation/Randomize.hpp
  QualityImprover/Relaxation/RelaxationSmoother.hpp
  QualityImprover/Relaxation/LaplacianSmoother.hpp
  QualityImprover/Relaxation/SmartLaplacianSmoother.hpp
  )

APPEND_SET(SOURCES
  QualityImprover/Relaxation/Randomize.cpp
  QualityImprover/Relaxation/RelaxationSmoother.cpp
  QualityImprover/Relaxation/LaplacianSmoother.cpp
  QualityImprover/Relaxation/SmartLaplacianSmoother.cpp
  )

#
# src/QualityMetric
#

APPEND_SET(HEADERS
  QualityMetric/AddQualityMetric.hpp
  QualityMetric/AveragingQM.hpp
  QualityMetric/EdgeQM.hpp
  QualityMetric/ElemSampleQM.hpp
  QualityMetric/ElementMaxQM.hpp
  QualityMetric/ElementAvgQM.hpp
  QualityMetric/ElementPMeanP.hpp
  QualityMetric/ElementQM.hpp
  QualityMetric/MultiplyQualityMetric.hpp
  QualityMetric/Debug/NumericalQM.hpp
  QualityMetric/Debug/CompareQM.hpp
  QualityMetric/PMeanPMetric.hpp
  QualityMetric/PowerQualityMetric.hpp
  QualityMetric/QualityMetric.hpp
  QualityMetric/ScalarAddQualityMetric.hpp
  QualityMetric/ScalarMultiplyQualityMetric.hpp
  QualityMetric/VertexMaxQM.hpp
  QualityMetric/VertexPMeanP.hpp
  QualityMetric/VertexQM.hpp 
  )

APPEND_SET(SOURCES
  QualityMetric/AddQualityMetric.cpp
  QualityMetric/AveragingQM.cpp
  QualityMetric/EdgeQM.cpp
  QualityMetric/ElemSampleQM.cpp
  QualityMetric/ElementMaxQM.cpp
  QualityMetric/ElementAvgQM.cpp
  QualityMetric/ElementPMeanP.cpp
  QualityMetric/ElementQM.cpp
  QualityMetric/MultiplyQualityMetric.cpp
  QualityMetric/Debug/NumericalQM.cpp
  QualityMetric/Debug/CompareQM.cpp
  QualityMetric/PMeanPMetric.cpp
  QualityMetric/PowerQualityMetric.cpp
  QualityMetric/QualityMetric.cpp
  QualityMetric/ScalarAddQualityMetric.cpp
  QualityMetric/ScalarMultiplyQualityMetric.cpp
  QualityMetric/VertexMaxQM.cpp
  QualityMetric/VertexPMeanP.cpp
  QualityMetric/VertexQM.cpp 
  )

#
# src/QualityMetric/Shape
#

APPEND_SET(HEADERS
  QualityMetric/Shape/AspectRatioGammaQualityMetric.hpp
  QualityMetric/Shape/ConditionNumberFunctions.hpp
  QualityMetric/Shape/ConditionNumberQualityMetric.hpp
  QualityMetric/Shape/IdealWeightInverseMeanRatio.hpp
  QualityMetric/Shape/IdealWeightMeanRatio.hpp
  QualityMetric/Shape/MeanRatioFunctions.hpp
  QualityMetric/Shape/VertexConditionNumberQualityMetric.hpp
  )

APPEND_SET(SOURCES
  QualityMetric/Shape/AspectRatioGammaQualityMetric.cpp
  QualityMetric/Shape/ConditionNumberQualityMetric.cpp
  QualityMetric/Shape/IdealWeightInverseMeanRatio.cpp
  QualityMetric/Shape/IdealWeightMeanRatio.cpp
  QualityMetric/Shape/VertexConditionNumberQualityMetric.cpp
  )

#
# src/QualityMetric/Smoothness
#

APPEND_SET(HEADERS
  QualityMetric/Smoothness/EdgeLengthQualityMetric.hpp
  QualityMetric/Smoothness/EdgeLengthRangeQualityMetric.hpp
  )

APPEND_SET(SOURCES
  QualityMetric/Smoothness/EdgeLengthQualityMetric.cpp
  QualityMetric/Smoothness/EdgeLengthRangeQualityMetric.cpp
  )

#
# src/QualityMetric/TMP
#

APPEND_SET(HEADERS
  QualityMetric/TMP/AffineMapMetric.hpp
  QualityMetric/TMP/TMPQualityMetric.hpp
  QualityMetric/TMP/AWQualityMetric.hpp 
  QualityMetric/TMP/TQualityMetric.hpp 
  )

APPEND_SET(SOURCES
  QualityMetric/TMP/AffineMapMetric.cpp
  QualityMetric/TMP/TargetMetricUtil.cpp
  QualityMetric/TMP/TargetMetricUtil.hpp
  QualityMetric/TMP/TMPQualityMetric.cpp 
  QualityMetric/TMP/AWQualityMetric.cpp 
  QualityMetric/TMP/TQualityMetric.cpp 
  )

#
# src/QualityMetric/Untangle
#

APPEND_SET(HEADERS
  QualityMetric/Untangle/UntangleBetaQualityMetric.hpp
  )

APPEND_SET(SOURCES
  QualityMetric/Untangle/UntangleBetaQualityMetric.cpp
  )

#
# src/QualityMetric/Volume
#

APPEND_SET(HEADERS
  QualityMetric/Volume/EdgeLengthMetric.hpp
  QualityMetric/Volume/LocalSizeQualityMetric.hpp
  QualityMetric/Volume/SizeMetric.hpp
  )

APPEND_SET(SOURCES
  QualityMetric/Volume/EdgeLengthMetric.cpp
  QualityMetric/Volume/LocalSizeQualityMetric.cpp
  QualityMetric/Volume/SizeMetric.cpp
  )


#
# src/TargetMetric
#

APPEND_SET(HEADERS
  TargetMetric/AWMetric.hpp
  TargetMetric/TMetric.hpp
  TargetMetric/TMetricBarrier.hpp
  TargetMetric/TMetricNonBarrier.hpp
  TargetMetric/AWMetricBarrier.hpp
  TargetMetric/AWMetricNonBarrier.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/TMPCommon.hpp
  TargetMetric/TMPDerivs.hpp
  TargetMetric/AWMetric.cpp
  TargetMetric/TMetric.cpp
  TargetMetric/TMetricBarrier.cpp
  TargetMetric/TMetricNonBarrier.cpp
  TargetMetric/AWMetricBarrier.cpp
  TargetMetric/AWMetricNonBarrier.cpp
  )


#
# src/TargetMetric/Misc
#

APPEND_SET(HEADERS
  TargetMetric/Misc/InvTransBarrier.hpp
  TargetMetric/Misc/TMixed.hpp
  TargetMetric/Misc/TOffset.hpp
  TargetMetric/Misc/TPower2.hpp
  TargetMetric/Misc/TScale.hpp
  TargetMetric/Misc/TSquared.hpp 
  TargetMetric/Misc/TSum.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/Misc/InvTransBarrier.cpp
  TargetMetric/Misc/TMixed.cpp
  TargetMetric/Misc/TOffset.cpp
  TargetMetric/Misc/TPower2.cpp
  TargetMetric/Misc/TScale.cpp
  TargetMetric/Misc/TSquared.cpp 
  TargetMetric/Misc/TSum.cpp
  )

#
# src/TargetMetric/Shape
#

APPEND_SET(HEADERS
  TargetMetric/Shape/AWShape2DB1.hpp
  TargetMetric/Shape/AWShape2DNB1.hpp
  TargetMetric/Shape/AWShape2DNB2.hpp
  TargetMetric/Shape/TInverseMeanRatio.hpp
  TargetMetric/Shape/TShape2DNB2.hpp
  TargetMetric/Shape/TShape3DB2.hpp
  TargetMetric/Shape/TShapeB1.hpp
  TargetMetric/Shape/TShapeNB1.hpp 
  )

APPEND_SET(SOURCES
  TargetMetric/Shape/AWShape2DB1.cpp
  TargetMetric/Shape/AWShape2DNB1.cpp
  TargetMetric/Shape/AWShape2DNB2.cpp
  TargetMetric/Shape/TInverseMeanRatio.cpp
  TargetMetric/Shape/TShape2DNB2.cpp
  TargetMetric/Shape/TShape3DB2.cpp
  TargetMetric/Shape/TShapeB1.cpp
  TargetMetric/Shape/TShapeNB1.cpp 
  )

#
# src/TargetMetric/ShapeOrient
#

APPEND_SET(HEADERS
  TargetMetric/ShapeOrient/AWShapeOrientNB1.hpp 
  TargetMetric/ShapeOrient/TShapeOrientB1.hpp 
  TargetMetric/ShapeOrient/TShapeOrientB2.hpp 
  TargetMetric/ShapeOrient/TShapeOrientNB1.hpp
  TargetMetric/ShapeOrient/TShapeOrientNB2.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/ShapeOrient/AWShapeOrientNB1.cpp 
  TargetMetric/ShapeOrient/TShapeOrientB1.cpp 
  TargetMetric/ShapeOrient/TShapeOrientB2.cpp 
  TargetMetric/ShapeOrient/TShapeOrientNB1.cpp
  TargetMetric/ShapeOrient/TShapeOrientNB2.cpp
  )

#
# src/TargetMetric/ShapeSize
#

APPEND_SET(HEADERS
  TargetMetric/ShapeSize/AWShapeSizeB1.hpp
  TargetMetric/ShapeSize/TShapeSize2DB2.hpp
  TargetMetric/ShapeSize/TShapeSize2DNB1.hpp
  TargetMetric/ShapeSize/TShapeSize2DNB2.hpp
  TargetMetric/ShapeSize/TShapeSize3DB2.hpp
  TargetMetric/ShapeSize/TShapeSize3DB4.hpp
  TargetMetric/ShapeSize/TShapeSize3DNB1.hpp
  TargetMetric/ShapeSize/TShapeSizeB1.hpp
  TargetMetric/ShapeSize/TShapeSizeB3.hpp
  TargetMetric/ShapeSize/TShapeSizeNB3.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/ShapeSize/AWShapeSizeB1.cpp
  TargetMetric/ShapeSize/TShapeSize2DB2.cpp
  TargetMetric/ShapeSize/TShapeSize2DNB1.cpp
  TargetMetric/ShapeSize/TShapeSize2DNB2.cpp
  TargetMetric/ShapeSize/TShapeSize3DB2.cpp
  TargetMetric/ShapeSize/TShapeSize3DB4.cpp
  TargetMetric/ShapeSize/TShapeSize3DNB1.cpp
  TargetMetric/ShapeSize/TShapeSizeB1.cpp
  TargetMetric/ShapeSize/TShapeSizeB3.cpp
  TargetMetric/ShapeSize/TShapeSizeNB3.cpp
  )

#
# src/TargetMetric/ShapeSizeOrient
#

APPEND_SET(HEADERS
  TargetMetric/ShapeSizeOrient/AWShapeSizeOrientNB1.hpp
  TargetMetric/ShapeSizeOrient/TShapeSizeOrientB1.hpp
  TargetMetric/ShapeSizeOrient/TShapeSizeOrientB2.hpp
  TargetMetric/ShapeSizeOrient/TShapeSizeOrientNB1.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/ShapeSizeOrient/AWShapeSizeOrientNB1.cpp
  TargetMetric/ShapeSizeOrient/TShapeSizeOrientB1.cpp
  TargetMetric/ShapeSizeOrient/TShapeSizeOrientB2.cpp
  TargetMetric/ShapeSizeOrient/TShapeSizeOrientNB1.cpp
  )

#
# src/TargetMetric/Size
#

APPEND_SET(HEADERS
  TargetMetric/Size/AWSizeB1.hpp
  TargetMetric/Size/AWSizeNB1.hpp
  TargetMetric/Size/TSizeB1.hpp
  TargetMetric/Size/TSizeNB1.hpp
  TargetMetric/Size/TTau.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/Size/AWSizeB1.cpp
  TargetMetric/Size/AWSizeNB1.cpp
  TargetMetric/Size/TSizeB1.cpp
  TargetMetric/Size/TSizeNB1.cpp
  TargetMetric/Size/TTau.cpp
  )

#
# src/TargetMetric/Untangle
#

APPEND_SET(HEADERS
  TargetMetric/Untangle/AWUntangleBeta.hpp
  TargetMetric/Untangle/TUntangle1.hpp
  TargetMetric/Untangle/TUntangleBeta.hpp
  TargetMetric/Untangle/TUntangleMu.hpp
  )

APPEND_SET(SOURCES
  TargetMetric/Untangle/AWUntangleBeta.cpp
  TargetMetric/Untangle/TUntangle1.cpp
  TargetMetric/Untangle/TUntangleBeta.cpp
  TargetMetric/Untangle/TUntangleMu.cpp
  )

#
# src/TargetCalculator
#

APPEND_SET(HEADERS
  TargetCalculator/CachingTargetCalculator.hpp
  TargetCalculator/IdealShapeTarget.hpp
  TargetCalculator/InverseMetricWeight.hpp
  TargetCalculator/JacobianCalculator.hpp
  TargetCalculator/LambdaTarget.hpp
  TargetCalculator/LambdaConstant.hpp
  TargetCalculator/LVQDTargetCalculator.hpp
  TargetCalculator/MetricWeight.hpp
  TargetCalculator/RefMeshTargetCalculator.hpp
  TargetCalculator/ReferenceMesh.hpp
  TargetCalculator/RefSizeTargetCalculator.hpp
  TargetCalculator/RemainingWeight.hpp
  TargetCalculator/TargetCalculator.hpp
  TargetCalculator/TargetReader.hpp
  TargetCalculator/TargetWriter.hpp
  TargetCalculator/TetDihedralWeight.hpp
  TargetCalculator/WeightCalculator.hpp
  TargetCalculator/WeightReader.hpp 
  )

APPEND_SET(SOURCES
  TargetCalculator/CachingTargetCalculator.cpp
  TargetCalculator/IdealShapeTarget.cpp
  TargetCalculator/InverseMetricWeight.cpp
  TargetCalculator/JacobianCalculator.cpp
  TargetCalculator/LambdaTarget.cpp
  TargetCalculator/LambdaConstant.cpp
  TargetCalculator/LVQDTargetCalculator.cpp
  TargetCalculator/MetricWeight.cpp
  TargetCalculator/RefMeshTargetCalculator.cpp
  TargetCalculator/ReferenceMesh.cpp
  TargetCalculator/RefSizeTargetCalculator.cpp
  TargetCalculator/RemainingWeight.cpp
  TargetCalculator/TargetCalculator.cpp
  TargetCalculator/TargetReader.cpp
  TargetCalculator/TargetWriter.cpp
  TargetCalculator/TetDihedralWeight.cpp
  TargetCalculator/WeightCalculator.cpp
  TargetCalculator/WeightReader.cpp 
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

