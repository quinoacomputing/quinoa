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


bool SubmaximalDF()
{
  int nx = 2;
  VectorType<double> vecType = new EpetraVectorType();

  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, 
    0.0, 1.0, nx, meshType);
  Mesh mesh = mesher.getMesh();

  CellFilter interior = new MaximalCellFilter();
  CellFilter pts = new DimensionalCellFilter(0);
  CellFilter middle = pts.subset(new CoordinateValueCellPredicate(0, 0.5));
    
  BasisFamily basis = new Lagrange(1);
  Expr u = new UnknownFunction(basis, "w");
  Expr v = new TestFunction(basis, "v");

  Expr grad = gradient(2);

  QuadratureFamily quad2 = new GaussianQuadrature(2);

  DiscreteSpace dsRes(mesh, basis, middle, vecType);
  Out::root() << "here's the map" << endl;
  dsRes.map()->print(Out::os());

  Expr uRes = new DiscreteFunction(dsRes, 1.0, "uRes");

  WatchFlag watch("watch");
  watch.setParam("discrete function evaluation", 5);
  watch.setParam("dof map setup", 5);
  watch.setParam("dof map access", 5);
  watch.setParam("eval mediator", 5);
  watch.setParam("integration", 1);
  watch.setParam("evaluation", 1);
  watch.deactivate();
  double f1 = L2Norm(mesh, middle, uRes, quad2, watch);

  Out::root() << "f1 = " << f1 << endl;

  bool test1 = std::fabs(f1 - std::sqrt(nx+1)) < 1.0e-10;
  Out::root() << "test1 pass=" << test1 << endl;

  Expr eqn = Integral(interior, (grad*u)*(grad*v), quad2);
  Expr bc = EssentialBC(middle, v*uRes*(u-uRes), quad2);

  LinearSolver<double> solver = LinearSolverBuilder::createSolver("amesos.xml");
    
  LinearProblem prob(mesh, eqn, bc, v, u, vecType);
  Expr u0 = prob.solve(solver);
  double f2 = L2Norm(mesh, interior, u0-1.0, quad2);
  Out::root() << "f2 = " << f2 << endl;

  bool test2 = std::fabs(f2) < 1.0e-10;
  Out::root() << "test2 pass=" << test2 << endl;

#ifdef CHECK_PROJECTION
  Expr y = new CoordExpr(1);
  Out::root() << "projecting by hand" << endl;

  watch.activate();
  Expr eqn_p = Integral(middle, v*(u-y), quad2, watch);
  Expr bc_p;

  LinearProblem prob_p(mesh, eqn_p, bc_p, v, u, vecType);
  Expr y_p = prob_p.solve(solver);



  Out::root() << "projecting by L2Proj obj" << endl;
  L2Projector proj(dsRes, y);
  Expr yRes = proj.project();

  double f3 = L2Norm(mesh, middle, yRes, quad2, watch);
  Out::root() << "f3 = " << f3 << endl;
  double g3 = 0.0;
  for (int i=0; i<=nx; i++) 
  {
    double z = i/(nx+1.0);
    g3 += z*z;
  }
  g3 = std::sqrt(g3);

  bool test3 = std::fabs(f3-g3) < 1.0e-10;
  Out::root() << "test3 pass=" << test3 << endl;
#endif

  return SundanceGlobal::checkTest(!test1 || !test2, 0.1);
}

