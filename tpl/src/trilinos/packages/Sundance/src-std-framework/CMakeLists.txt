

#
# A) Package-specific configuration options
#



#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})


#
# Elements files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../Playa/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-core/UserExprs)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-core/InternalExprs)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-core/FrameworkInterface)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Elements)

APPEND_SET(HEADERS
  Elements/SundanceBasisDOFTopologyBase.hpp
  Elements/SundanceBasisFamilyBase.hpp
  Elements/SundanceBasisFamily.hpp
  Elements/SundanceBasisReferenceEvaluationBase.hpp
  Elements/SundanceBernstein.hpp
  Elements/SundanceBubble.hpp
  Elements/SundanceClosedNewtonCotes.hpp
  Elements/SundanceClosedNewtonCotesType.hpp
  Elements/SundanceCubicHermite.hpp
  Elements/SundanceDiscreteFunction.hpp
  Elements/SundanceDiscreteFunctionData.hpp
  Elements/SundanceDiscreteSpace.hpp
  Elements/SundanceEdgeLocalizedBasis.hpp
  Elements/SundanceFeketeBrickQuadrature.hpp
  Elements/SundanceFeketeQuadrature.hpp
  Elements/SundanceFeketeQuadQuadrature.hpp
  Elements/SundanceFeketeQuadratureType.hpp
  Elements/SundanceFeketeTriangleQuadrature.hpp
  Elements/SundanceFourier.hpp
  Elements/SundanceTrapesoidQuadrature.hpp
  Elements/SundanceGaussLobattoQuadrature.hpp
  Elements/SundancePolygonQuadrature.hpp
  Elements/SundanceFuncWithBasis.hpp
  Elements/SundanceGauss1D.hpp
  Elements/SundanceGaussLobatto1D.hpp
  Elements/SundanceGaussianQuadrature.hpp
  Elements/SundanceGaussianQuadratureType.hpp
  Elements/SundanceLagrange.hpp
  Elements/SundanceNedelec.hpp
  Elements/SundanceQuadratureFamily.hpp
  Elements/SundanceQuadratureFamilyBase.hpp
  Elements/SundanceQuadratureType.hpp
  Elements/SundanceQuadratureTypeBase.hpp
  Elements/SundanceRaviartThomas.hpp
  Elements/SundanceReducedQuadrature.hpp
  Elements/SundanceTensorBasisBase.hpp
  Elements/SundanceTestFunction.hpp
  Elements/SundanceTestFunctionData.hpp
  Elements/SundanceTetQuadrature.hpp
  Elements/SundanceTriangleQuadrature.hpp
  Elements/SundanceQuadQuadrature.hpp
  Elements/SundanceBrickQuadrature.hpp
  Elements/SundanceUnknownFunction.hpp
  Elements/SundanceUnknownFunctionData.hpp
  Elements/SundanceLegendre.hpp
  Elements/SundanceSurfQuadrature.hpp
  )

APPEND_SET(SOURCES
  Elements/SundanceBasisDOFTopologyBase.cpp
  Elements/SundanceBasisFamily.cpp
  Elements/SundanceBasisFamilyBase.cpp
  Elements/SundanceBernstein.cpp
  Elements/SundanceBubble.cpp
  Elements/SundanceClosedNewtonCotes.cpp
  Elements/SundanceClosedNewtonCotesType.cpp
  Elements/SundanceCubicHermite.cpp
  Elements/SundanceDiscreteFunction.cpp
  Elements/SundanceDiscreteFunctionData.cpp
  Elements/SundanceDiscreteSpace.cpp
  Elements/SundanceEdgeLocalizedBasis.cpp
  Elements/SundanceFeketeBrickQuadrature.cpp
  Elements/SundanceFeketeQuadQuadrature.cpp
  Elements/SundanceFeketeQuadrature.cpp
  Elements/SundanceFeketeQuadratureType.cpp
  Elements/SundanceFeketeTriangleQuadrature.cpp
  Elements/SundanceFourier.cpp
  Elements/SundanceTrapesoidQuadrature.cpp
  Elements/SundanceGaussLobattoQuadrature.cpp
  Elements/SundancePolygonQuadrature.cpp
  Elements/SundanceFuncWithBasis.cpp
  Elements/SundanceGauss1D.cpp
  Elements/SundanceGaussLobatto1D.cpp
  Elements/SundanceGaussianQuadrature.cpp
  Elements/SundanceGaussianQuadratureType.cpp
  Elements/SundanceLagrange.cpp
  Elements/SundanceNedelec.cpp
  Elements/SundanceQuadratureFamily.cpp
  Elements/SundanceQuadratureFamilyBase.cpp
  Elements/SundanceQuadratureType.cpp
  Elements/SundanceRaviartThomas.cpp
  Elements/SundanceReducedQuadrature.cpp
  Elements/SundanceTestFunction.cpp
  Elements/SundanceTestFunctionData.cpp
  Elements/SundanceTetQuadrature.cpp
  Elements/SundanceTriangleQuadrature.cpp
  Elements/SundanceQuadQuadrature.cpp
  Elements/SundanceBrickQuadrature.cpp
  Elements/SundanceUnknownFunction.cpp
  Elements/SundanceUnknownFunctionData.cpp
  Elements/SundanceLegendre.cpp
  Elements/SundanceSurfQuadrature.cpp
  )

#
# DOF Mapping files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/DOFMapping)

APPEND_SET(HEADERS
  DOFMapping/SundanceBlock.hpp
  DOFMapping/SundanceDOFMapBase.hpp
  DOFMapping/SundanceDOFMapBuilder.hpp
  DOFMapping/SundanceInhomogeneousNodalDOFMap.hpp
  DOFMapping/SundanceInhomogeneousEdgeLocalizedDOFMap.hpp
  DOFMapping/SundanceMapStructure.hpp
  DOFMapping/SundanceMixedDOFMap.hpp
  DOFMapping/SundanceMixedDOFMapHN.hpp
  DOFMapping/SundanceNodalDOFMap.hpp
  DOFMapping/SundanceNodalDOFMapHN.hpp
  DOFMapping/SundanceInhomogeneousDOFMapHN.hpp
  DOFMapping/SundancePartialElementDOFMap.hpp
  DOFMapping/SundanceSpatiallyHomogeneousDOFMapBase.hpp
  DOFMapping/SundanceSubmaximalNodalDOFMap.hpp
  DOFMapping/SundanceHNDoFMapBase.hpp
  DOFMapping/SundanceHNDoFMapBaseHomogeneous.hpp
  DOFMapping/SundanceNodalDOFMapHN.hpp
  DOFMapping/SundanceUnfoldPeriodicDF.hpp
  )

APPEND_SET(SOURCES
  DOFMapping/SundanceDOFMapBase.cpp
  DOFMapping/SundanceDOFMapBuilder.cpp
  DOFMapping/SundanceInhomogeneousNodalDOFMap.cpp
  DOFMapping/SundanceInhomogeneousEdgeLocalizedDOFMap.cpp
  DOFMapping/SundanceMapStructure.cpp
  DOFMapping/SundanceMixedDOFMap.cpp
  DOFMapping/SundanceMixedDOFMapHN.cpp
  DOFMapping/SundanceNodalDOFMap.cpp
  DOFMapping/SundanceNodalDOFMapHN.cpp
  DOFMapping/SundanceInhomogeneousDOFMapHN.cpp
  DOFMapping/SundancePartialElementDOFMap.cpp
  DOFMapping/SundanceSpatiallyHomogeneousDOFMapBase.cpp
  DOFMapping/SundanceSubmaximalNodalDOFMap.cpp
  DOFMapping/SundanceUnfoldPeriodicDF.cpp
  )

#
# Assembly files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Assembly)
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Assembly/ElementIntegration)

APPEND_SET(HEADERS
  Assembly/SundanceAssembler.hpp
  Assembly/SundanceAssemblyKernelBase.hpp
  Assembly/SundanceElementIntegral.hpp
  Assembly/SundanceFunctionalAssemblyKernel.hpp
  Assembly/SundanceFunctionalGradientAssemblyKernel.hpp
  Assembly/SundanceFunctionalEvaluator.hpp
  Assembly/SundanceGrouperBase.hpp
  Assembly/SundanceIntegralGroup.hpp
  Assembly/SundanceIntegrationCellSpecifier.hpp
  Assembly/SundanceLocalDOFMap.hpp
  Assembly/SundanceLocalMatrixContainer.hpp
  Assembly/SundanceMapBundle.hpp
  Assembly/SundanceMatrixVectorAssemblyKernel.hpp
  Assembly/SundanceMaximalQuadratureIntegral.hpp
  Assembly/SundanceQuadratureEvalMediator.hpp
  Assembly/SundanceQuadratureIntegral.hpp
  Assembly/SundanceQuadratureIntegralBase.hpp
  Assembly/SundanceReducedIntegral.hpp
  Assembly/SundanceRefIntegral.hpp
  Assembly/SundanceStdFwkEvalMediator.hpp
  Assembly/SundanceTrivialGrouper.hpp
  Assembly/SundanceVectorAssemblyKernel.hpp
  Assembly/SundanceVectorFillingAssemblyKernel.hpp
  Assembly/SundanceCurveIntegralCalc.hpp
  Assembly/SundanceCurveQuadratureIntegral.hpp
  Assembly/SundanceCurveEvalMediator.hpp
  Assembly/SundanceSurf3DCalc.hpp
  )


#  Assembly/ElementIntegration/SundanceElementIntegralBase.hpp
#  Assembly/ElementIntegration/SundanceElementIntegralFunctional.hpp  
#  Assembly/ElementIntegration/SundanceElementIntegralFunctionalCell.hpp
#  Assembly/ElementIntegration/SundanceElementIntegralFunctionalFacet.hpp
#  Assembly/ElementIntegration/SundanceElementIntegralLinearForm.hpp
#  Assembly/ElementIntegration/SundanceElementIntegralLinearFormCell.hpp
#  Assembly/ElementIntegration/SundanceElementIntegralLinearFormFacet.hpp
#  Assembly/ElementIntegration/SundanceElementIntegralBilinearForm.hpp
#  Assembly/ElementIntegration/SundanceElementIntegralBilinearFormCell.hpp
#  Assembly/ElementIntegration/SundanceElementIntegralBilinearFormFacet.hpp
#  Assembly/ElementIntegration/SundanceIQI_HdivLF_DivV_Cell.hpp
#  Assembly/ElementIntegration/SundanceIQI_HdivLF_VN_Facet.hpp
#  Assembly/ElementIntegration/SundanceIQI_VecLF_Vi_Cell.hpp
#  Assembly/ElementIntegration/SundanceIQI_VecLF_Vi_Facet.hpp
#  Assembly/ElementIntegration/SundanceIQI_VecVecBLF_UiVj_Cell.hpp
#  Assembly/ElementIntegration/SundanceIQI_VecVecBLF_UiVj_Facet.hpp
#  Assembly/ElementIntegration/SundanceIQI_ScalVecBLF_UVj_Cell.hpp
#  Assembly/ElementIntegration/SundanceIQI_ScalVecBLF_UVj_Facet.hpp
#  Assembly/ElementIntegration/SundanceIQI_HdivScalBLF_DivUV_Cell.hpp
#  Assembly/ElementIntegration/SundanceIQI_HdivScalBLF_UNV_Facet.hpp
#  Assembly/ElementIntegration/SundanceIQI_HdivHdivBLF_DivUDivV_Cell.hpp

APPEND_SET(SOURCES
  Assembly/SundanceAssembler.cpp
  Assembly/SundanceElementIntegral.cpp
  Assembly/SundanceFunctionalAssemblyKernel.cpp
  Assembly/SundanceFunctionalEvaluator.cpp
  Assembly/SundanceGrouperBase.cpp
  Assembly/SundanceIntegralGroup.cpp
  Assembly/SundanceLocalDOFMap.cpp
  Assembly/SundanceLocalMatrixContainer.cpp
  Assembly/SundanceMapBundle.cpp
  Assembly/SundanceMatrixVectorAssemblyKernel.cpp
  Assembly/SundanceMaximalQuadratureIntegral.cpp
  Assembly/SundanceQuadratureEvalMediator.cpp
  Assembly/SundanceQuadratureIntegral.cpp
  Assembly/SundanceQuadratureIntegralBase.cpp
  Assembly/SundanceReducedIntegral.cpp
  Assembly/SundanceRefIntegral.cpp
  Assembly/SundanceStdFwkEvalMediator.cpp
  Assembly/SundanceTrivialGrouper.cpp
  Assembly/SundanceVectorAssemblyKernel.cpp
  Assembly/SundanceVectorFillingAssemblyKernel.cpp
  Assembly/SundanceCurveIntegralCalc.cpp
  Assembly/SundanceCurveQuadratureIntegral.cpp
  Assembly/SundanceCurveEvalMediator.cpp
  Assembly/SundanceSurf3DCalc.cpp
  )
#  Assembly/ElementIntegration/SundanceIQI_HdivLF_DivV_Cell.cpp
#  Assembly/ElementIntegration/SundanceIQI_HdivLF_VN_Facet.cpp
#  Assembly/ElementIntegration/SundanceIQI_VecLF_Vi_Cell.cpp
#  Assembly/ElementIntegration/SundanceIQI_VecLF_Vi_Facet.cpp
#  Assembly/ElementIntegration/SundanceIQI_VecVecBLF_UiVj_Cell.cpp
#  Assembly/ElementIntegration/SundanceIQI_VecVecBLF_UiVj_Facet.cpp
#  Assembly/ElementIntegration/SundanceIQI_ScalVecBLF_UVj_Cell.cpp
#  Assembly/ElementIntegration/SundanceIQI_ScalVecBLF_UVj_Facet.cpp
#  Assembly/ElementIntegration/SundanceIQI_HdivScalBLF_DivUV_Cell.cpp
#  Assembly/ElementIntegration/SundanceIQI_HdivScalBLF_UNV_Facet.cpp
#  Assembly/ElementIntegration/SundanceIQI_HdivHdivBLF_DivUDivV_Cell.cpp

#
# Problem files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Problem)

APPEND_SET(HEADERS
  Problem/Sundance.hpp
  Problem/SundanceAToCDensitySampler.hpp
  Problem/SundanceAToCPointLocator.hpp
  Problem/SundanceCToAInterpolator.hpp
  Problem/SundanceCoordinateSystem.hpp
  Problem/SundanceDoublingStepController.hpp
  Problem/SundanceEigensolution.hpp
  Problem/SundanceEventDetector.hpp
  Problem/SundanceExprComparison.hpp
  Problem/SundanceExprFieldWrapper.hpp
  Problem/SundanceFunctional.hpp
  Problem/SundanceL2Projector.hpp
  Problem/SundanceLinearEigenproblem.hpp
  Problem/SundanceLinearProblem.hpp
  Problem/SundanceLinearSolveDriver.hpp
  Problem/SundanceMeshIOUtils.hpp
  Problem/SundanceNitscheBC.hpp
  Problem/SundanceNLOp.hpp
  Problem/SundancePCDPreconditioner.hpp
  Problem/SundanceProblemTesting.hpp
  Problem/SundanceRivaraDriver.hpp
  Problem/SundanceStochBlockJacobiSolver.hpp
  Problem/SundanceNonlinearProblem.hpp 
  Problem/SundanceTransientStepProblem.hpp 
  Problem/SundanceUniformRefinementPair.hpp 
  )

APPEND_SET(SOURCES
  Problem/Sundance.cpp
  Problem/SundanceAToCDensitySampler.cpp
  Problem/SundanceAToCPointLocator.cpp
  Problem/SundanceCToAInterpolator.cpp
  Problem/SundanceCoordinateSystem.cpp
  Problem/SundanceDoublingStepController.cpp
  Problem/SundanceEventDetector.cpp
  Problem/SundanceExprComparison.cpp
  Problem/SundanceExprFieldWrapper.cpp
  Problem/SundanceFunctional.cpp
  Problem/SundanceL2Projector.cpp
  Problem/SundanceLinearEigenproblem.cpp
  Problem/SundanceLinearProblem.cpp
  Problem/SundanceLinearSolveDriver.cpp
  Problem/SundanceMeshIOUtils.cpp
  Problem/SundanceNitscheBC.cpp
  Problem/SundanceNLOp.cpp
  Problem/SundancePCDPreconditioner.cpp
  Problem/SundanceProblemTesting.cpp
  Problem/SundanceRivaraDriver.cpp
  Problem/SundanceStochBlockJacobiSolver.cpp
  Problem/SundanceNonlinearProblem.cpp 
  Problem/SundanceTransientStepProblem.cpp  
  Problem/SundanceUniformRefinementPair.cpp 
  )


#
# Element Matrix & Vector Transformation files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ElementTransformation)

APPEND_SET(HEADERS
  ElementTransformation/SundanceTransformationBase.hpp
  ElementTransformation/SundanceAssemblyTransformationBuilder.hpp
  ElementTransformation/SundanceDiscreteSpaceTransfBuilder.hpp
  ElementTransformation/SundanceTransformationHN.hpp
  ElementTransformation/SundanceInequivalentElementTransformation.hpp
  ElementTransformation/SundanceMatrixStore.hpp
)

APPEND_SET(SOURCES
  ElementTransformation/SundanceTransformationBase.cpp
  ElementTransformation/SundanceAssemblyTransformationBuilder.cpp
  ElementTransformation/SundanceDiscreteSpaceTransfBuilder.cpp
  ElementTransformation/SundanceTransformationHN.cpp
  ElementTransformation/SundanceInequivalentElementTransformation.cpp
  ElementTransformation/SundanceMatrixStore.cpp  
  )


IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  sundanceStdFwk
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS sundanceCore sundanceStdMesh sundanceInterop sundanceUtils sundancePlaya
  )
