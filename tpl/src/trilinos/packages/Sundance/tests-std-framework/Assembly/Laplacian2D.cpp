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

#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "PlayaTabs.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceExpr.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaEpetraVectorType.hpp"

using namespace Playa;
using namespace Teuchos;
using namespace Sundance;

static Time& totalTimer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}


int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      int pMax = 2;
      int dim=2;

      verbosity<RefIntegral>() = 0;

      CellType cellType = TriangleCell;

      Point a = Point(0.0, 0.0);
      Point b = Point(1.0, 0.0);
      Point c = Point(0.0, 1.0);

//       Point a = Point(1.0, 1.0);
//       Point b = Point(1.2, 1.6);
//       Point c = Point(0.8, 1.3);
      CellJacobianBatch JBatch;
      JBatch.resize(1, dim, dim);
      double* J = JBatch.jVals(0);
      J[0] = b[0] - a[0];
      J[1] = c[0] - a[0];
      J[2] = b[1] - a[1];
      J[3] = c[1] - a[1];

      QuadratureFamily q4 = new GaussianQuadrature(4);
          
      for (int p=1; p<=pMax; p++)
        {
          std::cerr << std::endl << "---------- p = " << p << " --------------" << std::endl;
          Tabs tab;
          BasisFamily P = new Lagrange(p);
      
          RCP<Array<double> > A = rcp(new Array<double>());
          RCP<Array<double> > Bxx = rcp(new Array<double>());
          RCP<Array<double> > Byy = rcp(new Array<double>());

          Array<double> constCoeff = tuple(1.0, 1.0);

          Array<int> alpha = tuple(0,1);
          Array<int> beta = tuple(0,1);
          ParametrizedCurve curve = new DummyParametrizedCurve();
           MeshType meshType = new BasicSimplicialMeshType();
           MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);
           Mesh mesh = mesher.getMesh();
           QuadratureFamily quad_1 = new GaussianQuadrature(2);
           RCP<Array<int> > cellLIDs;

          RefIntegral ref(dim, dim, cellType, P, alpha, 1, P, beta,quad_1, 1, isInternalBdry , curve , mesh);
          QuadratureIntegral qxx(dim, dim, cellType, P, tuple(0), 1, 
                                 P, tuple(0), 1, q4, curve , mesh,isInternalBdry);
          QuadratureIntegral qyy(dim, dim, cellType, P, tuple(1), 1, 
                                 P, tuple(1), 1, q4, curve , mesh,isInternalBdry);

          int nq = qxx.nQuad();
          Array<double> varCoeff(nq, 1.0);

          std::cerr << "============================= Doing reference integration =========== " << std::endl;

          ref.transformTwoForm(JBatch, constCoeff, cellLIDs, A);
          std::cerr << "============================= Doing quad integration xx =========== " << std::endl;
          qxx.transformTwoForm(JBatch, &(varCoeff[0]), cellLIDs , Bxx);
          std::cerr << "============================= Doing quad integration yy =========== " << std::endl;
          qyy.transformTwoForm(JBatch, &(varCoeff[0]), cellLIDs , Byy);

          std::cerr << "============================= Done integration =========== " << std::endl;
          std::cerr << tab << "transformed reference element" << std::endl;
          std::cerr << tab << "{";
          for (int r=0; r<ref.nNodesTest(); r++)
            {
              if (r!=0) std::cerr << ", ";
              std::cerr << "{";
              for (int c=0; c<ref.nNodesUnk(); c++)
                {
                  if (c!=0) std::cerr << ", ";
                  std::cerr << (*A)[r + ref.nNodesTest()*c];
                }
              std::cerr << "}";
            }
          std::cerr << "}" << std::endl;

          std::cerr << tab << "transformed Q_xx" << std::endl;
          std::cerr << tab << "{";
          for (int r=0; r<qxx.nNodesTest(); r++)
            {
              if (r!=0) std::cerr << ", ";
              std::cerr << "{";
              for (int c=0; c<qxx.nNodesUnk(); c++)
                {
                  if (c!=0) std::cerr << ", ";
                  int i = r + qxx.nNodesTest()*c;
                  double lapl = (*Bxx)[i];
                  std::cerr << lapl;
                }
              std::cerr << "}";
            }
          std::cerr << "}" << std::endl;

          std::cerr << tab << "transformed Q_yy" << std::endl;
          std::cerr << tab << "{";
          for (int r=0; r<qxx.nNodesTest(); r++)
            {
              if (r!=0) std::cerr << ", ";
              std::cerr << "{";
              for (int c=0; c<qxx.nNodesUnk(); c++)
                {
                  if (c!=0) std::cerr << ", ";
                  int i = r + qxx.nNodesTest()*c;
                  double lapl = (*Byy)[i];
                  std::cerr << lapl;
                }
              std::cerr << "}";
            }
          std::cerr << "}" << std::endl;
        }

      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
}
