

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

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/../src-utils/SundanceConfig.h
  )

#
# UserExprs files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-utils/Utilities)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-std-mesh/Interface)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/UserExprs)

APPEND_SET(HEADERS
  UserExprs/SundanceCellDiameterExpr.hpp
  UserExprs/SundanceCurveNormExpr.hpp
  UserExprs/SundanceCellVectorExpr.hpp
  UserExprs/SundanceComplexExpr.hpp
  UserExprs/SundanceCoordExpr.hpp
  UserExprs/SundanceDerivative.hpp
  UserExprs/SundanceEssentialBC.hpp
  UserExprs/SundanceExpr.hpp
  UserExprs/SundanceFunctionalDerivative.hpp
  UserExprs/SundanceIntegral.hpp
  UserExprs/SundanceParameter.hpp
  UserExprs/SundancePointwiseUserDefFunctor.hpp
  UserExprs/SundanceStdMathOps.hpp
  UserExprs/SundanceUnknownParameter.hpp
  UserExprs/SundanceUserDefFunctor.hpp
  UserExprs/SundanceUserDefOp.hpp 
  UserExprs/SundanceVectorCalculus.hpp 
  UserExprs/SundanceCurveExpr.hpp
  UserExprs/SundanceFunctionalEvaluatorBase.hpp
  )

APPEND_SET(SOURCES
  UserExprs/SundanceCellDiameterExpr.cpp
  UserExprs/SundanceCurveNormExpr.cpp
  UserExprs/SundanceCellVectorExpr.cpp
  UserExprs/SundanceComplexExpr.cpp
  UserExprs/SundanceCoordExpr.cpp
  UserExprs/SundanceDerivative.cpp
  UserExprs/SundanceEssentialBC.cpp
  UserExprs/SundanceExpr.cpp
  UserExprs/SundanceFunctionalDerivative.cpp
  UserExprs/SundanceIntegral.cpp
  UserExprs/SundanceParameter.cpp
  UserExprs/SundancePointwiseUserDefFunctor.cpp
  UserExprs/SundanceStdMathOps.cpp
  UserExprs/SundanceUnknownParameter.cpp
  UserExprs/SundanceUserDefFunctor.cpp
  UserExprs/SundanceUserDefOp.cpp
  UserExprs/SundanceVectorCalculus.cpp 
  )

#
# InternalExprs files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/InternalExprs)

APPEND_SET(HEADERS
  InternalExprs/SundanceBinaryExpr.hpp
  InternalExprs/SundanceConstantExpr.hpp
  InternalExprs/SundanceDerivOfSymbFunc.hpp
  InternalExprs/SundanceDiffOp.hpp
  InternalExprs/SundanceDiscreteFuncElement.hpp
  InternalExprs/SundanceEvaluatableExpr.hpp
  InternalExprs/SundanceExplicitFunctionalDerivativeElement.hpp
  InternalExprs/SundanceExprBase.hpp
  InternalExprs/SundanceExprWithChildren.hpp
  InternalExprs/SundanceFuncElementBase.hpp
  InternalExprs/SundanceFuncSetAccumulator.hpp
  InternalExprs/SundanceFunctionalPolynomial.hpp
  InternalExprs/SundanceFunctionWithID.hpp
  InternalExprs/SundanceListExpr.hpp
  InternalExprs/SundanceNonlinearExpr.hpp
  InternalExprs/SundanceNonlinearUnaryOp.hpp
  InternalExprs/SundanceParameterData.hpp
  InternalExprs/SundanceProductExpr.hpp
  InternalExprs/SundanceScalarExpr.hpp
  InternalExprs/SundanceSparsitySubset.hpp
  InternalExprs/SundanceSparsitySuperset.hpp
  InternalExprs/SundanceSpatiallyConstantExpr.hpp
  InternalExprs/SundanceSumExpr.hpp
  InternalExprs/SundanceSumOfBCs.hpp
  InternalExprs/SundanceSumOfIntegrals.hpp
  InternalExprs/SundanceSymbPreprocessor.hpp
  InternalExprs/SundanceSymbolicFuncElement.hpp
  InternalExprs/SundanceSymbolicFunc.hpp
  InternalExprs/SundanceTestFuncElement.hpp
  InternalExprs/SundanceUnaryExpr.hpp
  InternalExprs/SundanceUnaryMinus.hpp
  InternalExprs/SundanceUnknownFuncElement.hpp
  InternalExprs/SundanceUnknownParameterElement.hpp
  InternalExprs/SundanceUserDefFunctorElement.hpp
  InternalExprs/SundanceUserDefOpElement.hpp
  InternalExprs/SundanceZeroExpr.hpp 
  )

APPEND_SET(SOURCES
  InternalExprs/SundanceBinaryExpr.cpp
  InternalExprs/SundanceConstantExpr.cpp
  InternalExprs/SundanceDerivOfSymbFunc.cpp
  InternalExprs/SundanceDiffOp.cpp
  InternalExprs/SundanceDiscreteFuncElement.cpp
  InternalExprs/SundanceEvaluatableExpr.cpp
  InternalExprs/SundanceExplicitFunctionalDerivativeElement.cpp
  InternalExprs/SundanceExprBase.cpp
  InternalExprs/SundanceExprWithChildren.cpp
  InternalExprs/SundanceFuncElementBase.cpp
  InternalExprs/SundanceFunctionalPolynomial.cpp
  InternalExprs/SundanceFunctionWithID.cpp
  InternalExprs/SundanceListExpr.cpp
  InternalExprs/SundanceNonlinearExpr.cpp
  InternalExprs/SundanceNonlinearUnaryOp.cpp
  InternalExprs/SundanceProductExpr.cpp
  InternalExprs/SundanceScalarExpr.cpp
  InternalExprs/SundanceSparsitySubset.cpp
  InternalExprs/SundanceSparsitySuperset.cpp
  InternalExprs/SundanceSpatiallyConstantExpr.cpp
  InternalExprs/SundanceSumExpr.cpp
  InternalExprs/SundanceSumOfBCs.cpp
  InternalExprs/SundanceSumOfIntegrals.cpp
  InternalExprs/SundanceSymbPreprocessor.cpp
  InternalExprs/SundanceSymbolicFuncElement.cpp
  InternalExprs/SundanceSymbolicFunc.cpp
  InternalExprs/SundanceTestFuncElement.cpp
  InternalExprs/SundanceUnaryExpr.cpp
  InternalExprs/SundanceUnaryMinus.cpp
  InternalExprs/SundanceUnknownFuncElement.cpp
  InternalExprs/SundanceUnknownParameterElement.cpp
  InternalExprs/SundanceUserDefFunctorElement.cpp
  InternalExprs/SundanceUserDefOpElement.cpp
  InternalExprs/SundanceZeroExpr.cpp
  )

#
# Evaluation files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Evaluation)

APPEND_SET(HEADERS
  Evaluation/SundanceBinaryEvaluator.hpp
  Evaluation/SundanceCellDiameterEvaluator.hpp
  Evaluation/SundanceCurveNormEvaluator.hpp
  Evaluation/SundanceCellVectorEvaluator.hpp
  Evaluation/SundanceChainRuleEvaluator.hpp
  Evaluation/SundanceChainRuleSum.hpp
  Evaluation/SundanceConstantEvaluator.hpp
  Evaluation/SundanceCoordExprEvaluator.hpp
  Evaluation/SundanceDerivOfSymbFuncEvaluator.hpp
  Evaluation/SundanceDiscreteFuncEvaluator.hpp
  Evaluation/SundanceDiffOpEvaluator.hpp
  Evaluation/SundanceEFDEEvaluator.hpp
  Evaluation/SundanceEvalManager.hpp
  Evaluation/SundanceEvaluator.hpp
  Evaluation/SundanceEvaluatorFactory.hpp
  Evaluation/SundanceNonlinearUnaryOpEvaluator.hpp
  Evaluation/SundanceNullEvaluator.hpp
  Evaluation/SundanceProductEvaluator.hpp
  Evaluation/SundanceStringEvalMediator.hpp
  Evaluation/SundanceSubtypeEvaluator.hpp
  Evaluation/SundanceSumEvaluator.hpp
  Evaluation/SundanceSymbolicFuncEvaluator.hpp
  Evaluation/SundanceTempStack.hpp
  Evaluation/SundanceUnaryEvaluator.hpp
  Evaluation/SundanceUnaryMinusEvaluator.hpp
  Evaluation/SundanceUserDefOpCommonEvaluator.hpp
  Evaluation/SundanceUserDefOpEvaluator.hpp
  )

APPEND_SET(SOURCES
  Evaluation/SundanceCellDiameterEvaluator.cpp
  Evaluation/SundanceCurveNormEvaluator.cpp
  Evaluation/SundanceCellVectorEvaluator.cpp
  Evaluation/SundanceChainRuleEvaluator.cpp
  Evaluation/SundanceChainRuleSum.cpp
  Evaluation/SundanceConstantEvaluator.cpp
  Evaluation/SundanceCoordExprEvaluator.cpp
  Evaluation/SundanceDerivOfSymbFuncEvaluator.cpp
  Evaluation/SundanceDiscreteFuncEvaluator.cpp
  Evaluation/SundanceDiffOpEvaluator.cpp
  Evaluation/SundanceEFDEEvaluator.cpp
  Evaluation/SundanceEvalManager.cpp
  Evaluation/SundanceEvaluator.cpp
  Evaluation/SundanceNonlinearUnaryOpEvaluator.cpp
  Evaluation/SundanceNullEvaluator.cpp
  Evaluation/SundanceProductEvaluator.cpp
  Evaluation/SundanceStringEvalMediator.cpp
  Evaluation/SundanceSumEvaluator.cpp
  Evaluation/SundanceSymbolicFuncEvaluator.cpp
  Evaluation/SundanceTempStack.cpp
  Evaluation/SundanceUnaryMinusEvaluator.cpp
  Evaluation/SundanceUserDefOpCommonEvaluator.cpp
  Evaluation/SundanceUserDefOpEvaluator.cpp
  )

#
# FrameworkInterface files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/FrameworkInterface)

APPEND_SET(HEADERS
  FrameworkInterface/SundanceAbstractEvalMediator.hpp
  FrameworkInterface/SundanceAlgebraSpecifier.hpp
  FrameworkInterface/SundanceCommonFuncDataStub.hpp
  FrameworkInterface/SundanceComputationType.hpp
  FrameworkInterface/SundanceDeriv.hpp
  FrameworkInterface/SundanceDerivSet.hpp
  FrameworkInterface/SundanceDiscreteFuncDataStub.hpp
  FrameworkInterface/SundanceDiscreteFunctionStub.hpp
  FrameworkInterface/SundanceEquationSet.hpp
  FrameworkInterface/SundanceEvalContext.hpp
  FrameworkInterface/SundanceEvalVector.hpp
  FrameworkInterface/SundanceFunctionIdentifier.hpp
  FrameworkInterface/SundanceFunctionSupportResolver.hpp
  FrameworkInterface/SundanceLoadableVector.hpp
  FrameworkInterface/SundanceMultiIndex.hpp
  FrameworkInterface/SundanceMultipleDeriv.hpp
  FrameworkInterface/SundanceQuadratureFamilyStub.hpp
  FrameworkInterface/SundanceRegionQuadCombo.hpp
  FrameworkInterface/SundanceSpatialDerivSpecifier.hpp
  FrameworkInterface/SundanceSymbolicFuncDescriptor.hpp
  FrameworkInterface/SundanceTestFuncDataStub.hpp
  FrameworkInterface/SundanceTestFunctionStub.hpp
  FrameworkInterface/SundanceUnknownFuncDataStub.hpp
  FrameworkInterface/SundanceUnknownFunctionStub.hpp
  )

APPEND_SET(SOURCES
  FrameworkInterface/SundanceAbstractEvalMediator.cpp
  FrameworkInterface/SundanceAlgebraSpecifier.cpp
  FrameworkInterface/SundanceDeriv.cpp
  FrameworkInterface/SundanceDerivSet.cpp
  FrameworkInterface/SundanceDiscreteFunctionStub.cpp
  FrameworkInterface/SundanceEquationSet.cpp
  FrameworkInterface/SundanceEvalVector.cpp
  FrameworkInterface/SundanceFunctionIdentifier.cpp
  FrameworkInterface/SundanceFunctionSupportResolver.cpp
  FrameworkInterface/SundanceLoadableVector.cpp
  FrameworkInterface/SundanceMultiIndex.cpp
  FrameworkInterface/SundanceMultipleDeriv.cpp
  FrameworkInterface/SundanceQuadratureFamilyStub.cpp
  FrameworkInterface/SundanceRegionQuadCombo.cpp
  FrameworkInterface/SundanceSpatialDerivSpecifier.cpp
  FrameworkInterface/SundanceTestFunctionStub.cpp
  FrameworkInterface/SundanceUnknownFunctionStub.cpp
  )

#
# Transformations files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SymbolicTransformations)

APPEND_SET(HEADERS
  SymbolicTransformations/SundanceSymbolicTransformation.hpp
  SymbolicTransformations/SundanceSumTransformation.hpp
  SymbolicTransformations/SundanceSumTransformationSequence.hpp
  SymbolicTransformations/SundanceStdSumTransformations.hpp
  SymbolicTransformations/SundanceProductTransformation.hpp
  SymbolicTransformations/SundanceProductTransformationSequence.hpp
  SymbolicTransformations/SundanceStdProductTransformations.hpp
  )

APPEND_SET(SOURCES
  SymbolicTransformations/SundanceSymbolicTransformation.cpp
  SymbolicTransformations/SundanceSumTransformation.cpp
  SymbolicTransformations/SundanceSumTransformationSequence.cpp
  SymbolicTransformations/SundanceStdSumTransformations.cpp
  SymbolicTransformations/SundanceProductTransformation.cpp
  SymbolicTransformations/SundanceProductTransformationSequence.cpp
  SymbolicTransformations/SundanceStdProductTransformations.cpp
  )

#
# Spectral files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Spectral)

APPEND_SET(HEADERS
  Spectral/chaos.h
  Spectral/VECMAT.h
  Spectral/cijk.h
  Spectral/SundanceHermiteSpectralBasis.hpp
  Spectral/SundanceSpectralBasis.hpp
  Spectral/SundanceSpectralBasisBase.hpp
  Spectral/SundanceSpectralExpr.hpp 
  Spectral/SundanceSpectralPreprocessor.hpp 
  Spectral/SundanceStokhosBasisWrapper.hpp 
  )

APPEND_SET(SOURCES
  Spectral/SundanceHermiteSpectralBasis.cpp
  Spectral/SundanceSpectralExpr.cpp
  Spectral/SundanceSpectralPreprocessor.cpp 
  Spectral/SundanceStokhosBasisWrapper.cpp 
  )

#
# Test Utilities files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/TestingUtilities)

APPEND_SET(HEADERS
  TestingUtilities/SundanceADBasis.hpp
  TestingUtilities/SundanceADCoord.hpp
  TestingUtilities/SundanceADDerivative.hpp
  TestingUtilities/SundanceADField.hpp
  TestingUtilities/SundanceEvaluationTester.hpp
  TestingUtilities/SundanceTestDiscreteFuncData.hpp
  TestingUtilities/SundanceTestDiscreteFunction.hpp
  TestingUtilities/SundanceTestEvalMediator.hpp
  TestingUtilities/SundanceTestUnknownFuncData.hpp
  TestingUtilities/SundanceTestUnknownFunction.hpp
  )

APPEND_SET(SOURCES
  TestingUtilities/SundanceADBasis.cpp
  TestingUtilities/SundanceADCoord.cpp
  TestingUtilities/SundanceADDerivative.cpp
  TestingUtilities/SundanceADField.cpp
  TestingUtilities/SundanceEvaluationTester.cpp
  TestingUtilities/SundanceTestDiscreteFunction.cpp
  TestingUtilities/SundanceTestEvalMediator.cpp
  TestingUtilities/SundanceTestUnknownFuncData.cpp
  TestingUtilities/SundanceTestUnknownFunction.cpp
  )

#
# Parametrized curves
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ParametrizedCurves)

APPEND_SET(HEADERS
  ParametrizedCurves/SundanceCurveBase.hpp
  ParametrizedCurves/SundanceCircle.hpp
  ParametrizedCurves/SundanceBox2D.hpp
  ParametrizedCurves/SundanceBox3D.hpp
  ParametrizedCurves/SundanceEllipse2D.hpp
  ParametrizedCurves/SundanceSphere.hpp
  ParametrizedCurves/SundanceLine2D.hpp
  ParametrizedCurves/SundancePlane3D.hpp
  ParametrizedCurves/SundanceParamCurveIntegral.hpp
  ParametrizedCurves/SundanceParametrizedCurve.hpp
  ParametrizedCurves/SundanceDummyParametrizedCurve.hpp
  ParametrizedCurves/SundancePolygon2D.hpp
  ParametrizedCurves/SundanceTriangleSurf3D.hpp
  ParametrizedCurves/SundanceCurveCollection.hpp
  )

APPEND_SET(SOURCES
  ParametrizedCurves/SundanceCircle.cpp
  ParametrizedCurves/SundanceBox2D.cpp
  ParametrizedCurves/SundanceBox3D.cpp
  ParametrizedCurves/SundanceEllipse2D.cpp
  ParametrizedCurves/SundanceSphere.cpp
  ParametrizedCurves/SundanceLine2D.cpp
  ParametrizedCurves/SundancePlane3D.cpp
  ParametrizedCurves/SundanceParamCurveIntegral.cpp
  ParametrizedCurves/SundanceParametrizedCurve.cpp
  ParametrizedCurves/SundancePolygon2D.cpp
  ParametrizedCurves/SundanceTriangleSurf3D.cpp
  ParametrizedCurves/SundanceCurveCollection.cpp

  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  sundanceCore
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS sundanceUtils sundanceInterop sundancePlaya
  )
