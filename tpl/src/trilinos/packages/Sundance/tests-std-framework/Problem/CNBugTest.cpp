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


#include "Sundance.hpp"


// written by Michael Ulbrich, TU Muenchen, Dec 30, 2010

bool CNBugTest()
{
  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  MeshType meshType = new BasicSimplicialMeshType();
  int nx = 1;

  MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
  Mesh mesh = mesher.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter(); 

  /* Create unknown function */
  BasisFamily L1=new Lagrange(1);
  Expr u = new UnknownFunction(L1, "u");
  Expr w = new TestFunction(L1, "w");

  Expr dx = new Derivative(0);
      
  Expr x = new CoordExpr(0);

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad = new GaussianQuadrature(4);

  DiscreteSpace discSpaceL1(mesh, L1, vecType);
  Expr ud = x; //proj.project();
    
  WatchFlag watch("watch");
  watch.setParam("evaluation", 5);
  watch.setParam("evaluator setup", 5);
  watch.setParam("discrete function evaluation", 1);
  watch.setParam("integration setup", 0);
  watch.setParam("symbolic preprocessing", 3);
  watch.setParam("integration", 0);

  Expr eqn1 = Integral(interior, w*(dx*(u+ud)), quad, watch);
  Expr eqn2 = Integral(interior, w*(dx*u)+w*(dx*ud), 
    quad, watch);
  Expr bc;

  Out::root() << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" 
              << " @@@@@@@@@@@@@@@ creating BAD operator @@@@@@@@@@@@@\n" 
              << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" 
              << endl;
    
  LinearProblem prob1(mesh, eqn1, bc, w, u, vecType);
  LinearOperator<double> A1 = prob1.getOperator();

    
  Out::root() << "BAD operator = " << endl << A1 << endl << endl << endl;


  Out::root() << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" 
              << " @@@@@@@@@@@@@@@ creating GOOD operator @@@@@@@@@@@@\n" 
              << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
              << endl;
  LinearProblem prob2(mesh, eqn2, bc, w, u, vecType);
  LinearOperator<double> A2 = prob2.getOperator();
  Out::root() << "GOOD operator = " << endl << A2 << endl;

  Vector<double> b = A2.domain().createMember();
  b.randomize();

  Vector<double> r = A2*b - A1*b;

  Out::root() << "difference in operator application = " 
              << r.norm2() << endl;
    
  double tol = 1000.0;
  return SundanceGlobal::checkTest(r.norm2(), tol);

}
