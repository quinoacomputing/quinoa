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
  int stat = 0;
  int verb=1;
  try
  {
    GlobalMPISession session(&argc, &argv);

    TimeMonitor t(totalTimer());

    int pMax = 2;
    int dim=2;

    bool isInternalBdry = false;

    Utils::setChopVal(1.0e-14);

    CellType cellType = TriangleCell;

    //       Point a = Point(1.0, 1.0);
    //       Point b = Point(1.2, 1.6);
    //       Point c = Point(0.8, 1.3);

    Point a = Point(0.0, 0.0);
    Point b = Point(1.0, 0.0);
    Point c = Point(0.0, 1.0);

    Point d = Point(0.1, 0.1);
    Point e = Point(1.0, 0.0);
    Point f = Point(0.0, 1.0);

    int nCells = 2;

    CellJacobianBatch JBatch;
    JBatch.resize(nCells, 2, 2);
    double* J = JBatch.jVals(0);
    J[0] = b[0] - a[0];
    J[1] = c[0] - a[0];
    J[2] = b[1] - a[1];
    J[3] = c[1] - a[1];

    J[4] = e[0] - d[0];
    J[5] = f[0] - d[0];
    J[6] = e[1] - d[1];
    J[7] = f[1] - d[1];


      
    Array<int> dummy;
    double coeff = 1.0;
    RCP<Array<double> > A = rcp(new Array<double>());
    RCP<Array<double> > B = rcp(new Array<double>());

    QuadratureFamily q4 = new GaussianQuadrature(4);

    int nErrors = 0;

    std::cerr << std::endl << std::endl 
         << "---------------- One-forms --------------------" 
         << std::endl << std::endl;
    for (int p=0; p<=pMax; p++)
    {
      BasisFamily P = new Lagrange(p);
      for (int dp=0; dp<=1; dp++)
      {
        if (dp > p) continue;

        int numTestDir = 1;
        if (dp==1) numTestDir = dim;
        for (int t=0; t<numTestDir; t++)
        {
          int alpha = t;
          Tabs tab;

          ParametrizedCurve curve = new DummyParametrizedCurve();
          MeshType meshType = new BasicSimplicialMeshType();
          MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);
          Mesh mesh = mesher.getMesh();
          RCP<Array<int> > cellLIDs;

          RefIntegral ref(dim, cellType, dim, cellType, P, alpha, dp, q4 , isInternalBdry, curve, mesh ,verb);
          A->resize(JBatch.numCells() * ref.nNodes());
          for (int ai=0; ai<A->size(); ai++) (*A)[ai]=0.0;
          ref.transformOneForm(JBatch, JBatch, dummy, cellLIDs , coeff, A);
          std::cerr << tab << "transformed reference element" << std::endl;
          if (dp>0) std::cerr << tab << "test diff direction=" << t << std::endl;
          for (int cell=0; cell<nCells; cell++)
          {
            std::cerr << tab << "{";
            for (int r=0; r<ref.nNodesTest(); r++)
            {
              if (r!=0) std::cerr << ", ";
              std::cerr << Utils::chop((*A)[cell*ref.nNodesTest()+r]);
            }
            std::cerr << "}" << std::endl;
          }
          QuadratureIntegral quad(dim, cellType, dim, cellType, P, alpha, dp, q4, isInternalBdry, curve, mesh, verb);
          Array<double> quadCoeff(2*quad.nQuad(), 1.0);
          B->resize(JBatch.numCells() * quad.nNodes());
          for (int ai=0; ai<B->size(); ai++) (*B)[ai]=0.0;
          quad.transformOneForm(JBatch, JBatch, dummy, cellLIDs , &(quadCoeff[0]), B);
          std::cerr << tab << "transformed quad element" << std::endl;
          if (dp>0) std::cerr << tab << "test diff direction =" << t << std::endl;
          for (int cell=0; cell<nCells; cell++)
          {
            std::cerr << tab << "{";
            for (int r=0; r<quad.nNodesTest(); r++)
            {
              if (r!=0) std::cerr << ", ";
              std::cerr << Utils::chop((*B)[cell*ref.nNodesTest()+r]);
            }
            std::cerr << "}" << std::endl;
          }

          std::cerr << tab << "MISFIT quad-ref" << std::endl;
          std::cerr << tab << "test diff order =" << dp << std::endl;
          if (dp>0) std::cerr << tab << "test diff direction =" << t << std::endl;
          bool OK = true;
          for (int cell=0; cell<nCells; cell++)
          {
            std::cerr << tab << "{";
            for (int r=0; r<quad.nNodesTest(); r++)
            {
              if (r!=0) std::cerr << ", ";
              int i = cell*ref.nNodesTest()+r;
              double err = fabs(Utils::chop((*B)[i] - (*A)[i]));
              if (err > 1.0e-14) 
              {
                OK = false;
              }
              std::cerr << err;
            }
            std::cerr << "}" << std::endl;
          }
                  
          if (!OK) 
          {
            nErrors ++;
            std::cerr << "ERROR DETECTED!!! p=" << p
                 << "  t=" << t  << std::endl;
          }
        }
      }
    }
         




    std::cerr << std::endl << std::endl 
         << "---------------- Two-forms --------------------" 
         << std::endl << std::endl;

    for (int p=0; p<=pMax; p++)
    {
      BasisFamily P = new Lagrange(p);
      for (int dp=0; dp<=1; dp++)
      {
        if (dp > p) continue;
        int numTestDir = 1;
        if (dp==1) numTestDir = dim;
        for (int q=0; q<=pMax; q++)
        {
          BasisFamily Q = new Lagrange(q);
          for (int dq=0; dq<=1; dq++)
          {
            if (dq > q) continue;
            for (int t=0; t<numTestDir; t++)
            {
              int alpha = t;
              int numUnkDir = 1;
              if (dq==1) numUnkDir = dim;
              for (int u=0; u<numUnkDir; u++)
              {
                   ParametrizedCurve curve = new DummyParametrizedCurve();
                   MeshType meshType = new BasicSimplicialMeshType();
                   MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);
                   Mesh mesh = mesher.getMesh();
                   QuadratureFamily quad_1 = new GaussianQuadrature(2);
                   RCP<Array<int> > cellLIDs;

                Tabs tab;
                //                              if (p==0 || q==0 || dp==0 || dq==0 || u==1
                //  || t==1) continue;
                int beta = u;
                RefIntegral ref(dim, cellType, dim, cellType, P, alpha,
                  dp, Q, beta, dq, quad_1 , isInternalBdry, curve , mesh , verb);
                A->resize(JBatch.numCells() * ref.nNodes());
                for (int ai=0; ai<A->size(); ai++) (*A)[ai]=0.0;
                ref.transformTwoForm(JBatch, JBatch, dummy, cellLIDs , coeff, A);
                std::cerr << tab << "transformed ref element" << std::endl;
                std::cerr << tab << "test diff order = " << dp << std::endl;
                if (dp>0) std::cerr << tab << "t=dx(" << t << ")" << std::endl;
                std::cerr << tab << "unk diff order = " << dq << std::endl;
                if (dq>0) std::cerr << tab << "u=dx(" << u << ")" << std::endl;

                for (int cell=0; cell<nCells; cell++)
                {
                  std::cerr << tab << "cell=" << cell << " {";
                  for (int r=0; r<ref.nNodesTest(); r++)
                  {
                    if (r!=0) std::cerr << ", ";
                    std::cerr << "{";
                    for (int c=0; c<ref.nNodesUnk(); c++)
                    {
                      if (c!=0) std::cerr << ", ";
                      std::cerr << Utils::chop((*A)[r + ref.nNodesTest()*(c + cell*ref.nNodesUnk())]);
                    }
                    std::cerr << "}";
                  }
                  std::cerr << "}" << std::endl;
                }


                QuadratureIntegral quad(dim, cellType, dim, cellType, P, alpha,
                  dp, Q, beta, dq, q4, isInternalBdry,curve , mesh , verb);
                Array<double> quadCoeff(2*quad.nQuad(), 1.0);
                B->resize(JBatch.numCells() * quad.nNodes());
                for (int ai=0; ai<B->size(); ai++) (*B)[ai]=0.0;
                quad.transformTwoForm(JBatch, JBatch, dummy, cellLIDs , &(quadCoeff[0]), B);

                std::cerr << tab << "transformed quad element" << std::endl;
                std::cerr << tab << "test diff order = " << dp << std::endl;
                if (dp>0) std::cerr << tab << "t=dx(" << t << ")" << std::endl;
                std::cerr << tab << "unk diff order = " << dq << std::endl;
                if (dq>0) std::cerr << tab << "u=dx(" << u << ")" << std::endl;

                for (int cell=0; cell<nCells; cell++)
                {
                  std::cerr << tab << "cell=" << cell << " {";
                  for (int r=0; r<ref.nNodesTest(); r++)
                  {
                    if (r!=0) std::cerr << ", ";
                    std::cerr << "{";
                    for (int c=0; c<ref.nNodesUnk(); c++)
                    {
                      if (c!=0) std::cerr << ", ";
                      std::cerr << Utils::chop((*B)[r + ref.nNodesTest()*(c + cell*ref.nNodesUnk())]);
                    }
                    std::cerr << "}";
                  }
                  std::cerr << "}" << std::endl;   
                }

                bool OK = true;
                std::cerr << tab << "MISMATCH quad - ref" << std::endl;
                std::cerr << tab << "test diff order = " << dp << std::endl;
                if (dp>0) std::cerr << tab << "t=dx(" << t << ")" << std::endl;
                std::cerr << tab << "unk diff order = " << dq << std::endl;
                if (dq>0) std::cerr << tab << "u=dx(" << u << ")" << std::endl;

                for (int cell=0; cell<nCells; cell++)
                {
                  std::cerr << tab << "cell #" << cell << " {";
                              
                  for (int r=0; r<ref.nNodesTest(); r++)
                  {
                    if (r!=0) std::cerr << ", ";
                    std::cerr << "{";
                    for (int c=0; c<ref.nNodesUnk(); c++)
                    {
                      if (c!=0) std::cerr << ", ";
                      int i = r + ref.nNodesTest()*(c + cell*ref.nNodesUnk());
                      double err = fabs(Utils::chop((*B)[i] - (*A)[i]));
                      if (err > 1.0e-14) OK = false;
                      std::cerr << err;
                    }
                    std::cerr << "}";
                  }
                  std::cerr << "}" << std::endl;
                }
                if (!OK) 
                {
                  nErrors ++;
                  std::cerr << "ERROR DETECTED!!! p=" << p
                       << " dp=" << dp << "  t=" << t  
                       << " q=" << q << "  dq=" << dq
                       << "  u=" << u << std::endl;
                }

                std::cerr << std::endl << std::endl << std::endl << std::endl;
              }
            }
          }
        }
      }
    }

    std::cerr << "total quadrature flops: " << QuadratureIntegral::totalFlops() 
         << std::endl;
    std::cerr << "total ref integration flops: " << RefIntegral::totalFlops() 
         << std::endl;

    if (nErrors == 0)
    {
      std::cerr << "Transformed integral test PASSED" << std::endl;
    }
    else
    {
      stat = -1;
      std::cerr << "Transformed integral test FAILED" << std::endl;
    }
    TimeMonitor::summarize();
  }
	catch(std::exception& e)
  {
    stat = -1;
    std::cerr << "Transformed integral test FAILED" << std::endl;
    std::cerr << e.what() << std::endl;
  }

  return stat;
  
}
