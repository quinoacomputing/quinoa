/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCE_H
#define SUNDANCE_H

/* Utilities */
#include "SundanceDefs.hpp"
#include "SundanceOut.hpp"
#include "SundancePathUtils.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

/* Symbolics */
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceUnknownParameter.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceCurveNormExpr.hpp"
#include "SundanceCellVectorExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceParameter.hpp"
#include "SundancePointwiseUserDefFunctor.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceUserDefOp.hpp"
#include "SundanceVectorCalculus.hpp"
#include "SundanceCurveExpr.hpp"

/* Parametrized Curves */
#include "SundanceCircle.hpp"
#include "SundanceBox2D.hpp"
#include "SundanceBox3D.hpp"
#include "SundanceLine2D.hpp"
#include "SundancePlane3D.hpp"
#include "SundanceEllipse2D.hpp"
#include "SundanceSphere.hpp"
#include "SundancePolygon2D.hpp"
#include "SundanceTriangleSurf3D.hpp"
#include "SundanceDummyParametrizedCurve.hpp"
#include "SundanceParamCurveIntegral.hpp"
#include "SundanceParametrizedCurve.hpp"
#include "SundanceCurveIntegralCalc.hpp"
#include "SundanceCurveCollection.hpp"

/* Meshes */
#include "SundanceMesh.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceMeshTransformation.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceExtrusionMeshTransformation.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "SundanceSerialPartitionerBase.hpp"
#include "SundanceFileIOChacoPartitioner.hpp"
#include "SundanceTriangleMeshReader.hpp"
#include "SundanceExodusNetCDFMeshReader.hpp"
#include "SundanceExodusMeshReader.hpp"
#include "SundanceMeshBuilder.hpp"
#include "SundanceBamgMeshReader.hpp"
#include "SundanceHNMesher2D.hpp"
#include "SundanceHNMeshType2D.hpp"
#include "SundanceHNMesher3D.hpp"
#include "SundanceHNMeshType3D.hpp"


#ifdef HAVE_SUNDANCE_PEANO
	#ifdef HAVE_SUNDANCE_PEANO_NO_2D
	#else
		#include "SundancePeanoMeshType2D.hpp"
		#include "SundancePeanoMesher2D.hpp"
	#endif
	#ifdef HAVE_SUNDANCE_PEANO_NO_3D
	#else
		#include "SundancePeanoMeshType3D.hpp"
		#include "SundancePeanoMesher3D.hpp"
	#endif
#endif

/* Mesh refinement*/
#include "SundanceRefinementBase.hpp"
#include "SundanceRefinementClass.hpp"

/* Cell filters */
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceBoundaryCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceDomainDefinition.hpp"

/* Writers */
#include "SundanceFieldWriter.hpp"
#include "SundanceMatlabWriter.hpp" // deprecated
#include "SundanceDSVWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceTriangleWriter.hpp"
#include "SundanceVTKWriter.hpp"
#include "SundanceExodusWriter.hpp"


/* FE  */
#include "SundanceBasisFamily.hpp"
#include "SundanceBernstein.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceLegendre.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceFeketeQuadrature.hpp"
#include "SundanceTrapesoidQuadrature.hpp"
#include "SundanceGaussLobattoQuadrature.hpp"
#include "SundancePolygonQuadrature.hpp"
#include "SundanceSurfQuadrature.hpp"

/* Spectral */
#include "SundanceHermiteSpectralBasis.hpp"

/* Problem level classes */
#include "SundanceCoordinateSystem.hpp"
#include "SundanceLinearProblem.hpp"
#include "SundanceLinearEigenproblem.hpp"
#include "SundanceL2Projector.hpp"
#include "SundanceNonlinearProblem.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "SundanceFunctional.hpp"
#include "SundanceRivaraDriver.hpp"
#include "SundanceExprFieldWrapper.hpp"
#include "SundanceNitscheBC.hpp"
#include "SundanceBlock.hpp"


/* Solvers & stuff */
#include "PlayaVectorType.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaEpetraMatrixMatrixProduct.hpp"
#include "PlayaEpetraMatrixMatrixSum.hpp"
#include "PlayaEpetraMatrixOps.hpp"
#include "PlayaBICGSTABSolverDecl.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaMLOperator.hpp"
#include "PlayaParameterListPreconditionerFactory.hpp"
#include "SundancePCDPreconditioner.hpp"
#include "PlayaEpetraMatrixMatrixProduct.hpp"
#include "PlayaEpetraMatrixMatrixSum.hpp"
#include "PlayaEpetraMatrixOps.hpp"



/* Nonlinear solvers */
#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "PlayaNOXSolver.hpp"
#include "PlayaNonlinearSolver.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"


/* Eigensolvers */
#include "PlayaAnasaziEigensolverDecl.hpp"
#include "PlayaEigensolver.hpp"


/* Atomistic/continuum */
#include "SundanceAToCDensitySampler.hpp"
#include "SundanceCToAInterpolator.hpp"


/* do explicit qualification of List to avoid conflicts
 * with the unfriendly lack of namespaces in MPI2C++
 */

namespace Sundance
{

using namespace PlayaExprTemplates;
using namespace Playa;
using namespace Teuchos;
using Sundance::List;
  using std::setw;
  using std::exception;


/**
 * Class Sundance provides several static methods for
 * managing the environment of a simulation. Every simulation code
 * should begin with a call to Sundance::init() and end with
 * a call to Sundance::finalize().
 */
class SundanceGlobal
{
public:

  /** */
  static void setOption(const std::string& optionName,
    int& value,
    const std::string& helpMsg);

  /** */
  static void setOption(const std::string& optionName,
    std::string& value,
    const std::string& helpMsg);

  /** */
  static void setOption(const std::string& optionName,
    double& value,
    const std::string& helpMsg);

  /** */
  static void setOption(const std::string& optionTrueName,
    const std::string& optionFalseName,
    bool& value,
    const std::string& helpMsg);


  /** 
   * Do initialization steps such as starting MPI (if necessary), 
   * parsing the Unix command
   * line, and reading runtime options from the configuration file.
   * MPI is initialized through a call to Teuchos::GlobalMPISession, 
   * which in turn checks whether MPI needs initialization and calls
   * MPI_Init() if necessary. If some other library you're using has
   * its own MPI initialization system, it is thus perfectly safe to
   * call Sundance::init() as well.
   */
  static int init(int* argc, char*** argv);
    
  /** 
   * Do finalization steps such as calling MPI_finalize() (if necessary),
   * and collecting and printing timing information.
   */
  static int finalize();
    
  /** */
  static void handleException(std::exception& e);

  /** */
  static Teuchos::FancyOStream& os() ;


  /** */
  static bool passFailTest(bool pass);

  /** */
  static bool passFailTest(double error, double tol);

  /** */
  static bool passFailTest(const std::string& statusMsg,
    bool status, double error, double tol);


  /** Set to true if a message should be written by each processor
   * at startup. */
  static bool& showStartupMessage();

  /** Decide whether to skip timing outputs to work around
   * a trilinos 6.0.x bug */
  static bool& skipTimingOutput()
    {static bool rtn=false; return rtn;}

  /** */
  static int& testStatus() {static int rtn = -1; return rtn;}



  static CommandLineProcessor& clp()
    {static CommandLineProcessor rtn; return rtn;}

  static bool checkTest(double error, double tol);


};


/** \relates SundanceGlobal */
void handleException(std::exception& e);

/** \relates SundanceGlobal */
bool passFailTest(bool pass);

/** \relates SundanceGlobal */
bool passFailTest(double error, double tol);

/** \relates SundanceGlobal */
bool passFailTest(const std::string& statusMsg,
  bool status, double error, double tol);

/** \relates SundanceGlobal */
int& testStatus() ;

/** \relates SundanceGlobal */
CommandLineProcessor& clp() ;

/** \relates SundanceGlobal */
int init(int* argc, char*** argv);

/** \relates SundanceGlobal */
int finalize();


/** */
void setOption(const std::string& optionName,
  int& value,
  const std::string& helpMsg);

/** */
void setOption(const std::string& optionName,
  std::string& value,
  const std::string& helpMsg);

/** */
void setOption(const std::string& optionName,
  double& value,
  const std::string& helpMsg);

/** */
void setOption(const std::string& optionTrueName,
  const std::string& optionFalseName,
  bool& value,
  const std::string& helpMsg);

}




#endif

