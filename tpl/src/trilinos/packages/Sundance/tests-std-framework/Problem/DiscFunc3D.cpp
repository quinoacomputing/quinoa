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
#include "SundanceElementIntegral.hpp"
using Sundance::List;




/**
 * This program tests the ability to compute derivatives of discrete
 * functions in thre dimensions.
 */


#ifdef HAVE_SUNDANCE_EXODUS

bool DiscFunc3D()
{
  std::string meshFile="plate3D-0";
  int order = 2;

  VectorType<double> vecType = new EpetraVectorType();

  MeshType meshType = new BasicSimplicialMeshType();
      
  MeshSource meshSrc
    =  new ExodusMeshReader(meshFile, meshType);
  Mesh mesh = meshSrc.getMesh();

  CellFilter interior = new MaximalCellFilter();

    
  BasisFamily basis = new Lagrange(order);
  DiscreteSpace discSpace(mesh, basis, vecType);


  Expr dx = new Derivative(0);
  Expr dy = new Derivative(1);
  Expr dz = new Derivative(2);
  Expr grad = List(dx, dy, dz);

  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);
  Expr z = new CoordExpr(2);
    

  /* */
  Expr f = 0.5*(x*x + y*y + z*z);
  L2Projector proj(discSpace, f);
  Expr u0 = proj.project();
  L2Projector proj2(discSpace, dx*u0 - dx*f);
  Expr du0_dx = proj2.project();
  L2Projector proj3(discSpace, dy*u0 - dy*f);
  Expr du0_dy = proj3.project();
  L2Projector proj4(discSpace, dz*u0 - dz*f);
  Expr du0_dz = proj4.project();

  Expr F = List(x,y,z);
  Expr n = CellNormalExpr(3, "n");
  Expr v = new TestFunction(basis);
  Expr u = new UnknownFunction(basis);
  CellFilter bdry = new BoundaryCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(4);
    
  Expr eqn = Integral(interior, u*v  + F*(grad*v) + v*(grad*F), quad)
    - Integral(bdry, v*(F*n), quad);
  Expr bc;
  LinearProblem prob(mesh, eqn, bc, v, u, vecType);
    
  LinearSolver<double> solver 
    = LinearSolverBuilder::createSolver("amesos.xml");

  Expr soln = prob.solve(solver);

  FieldWriter w = new VTKWriter("DiscFunc3D");
  w.addMesh(mesh);
  w.addField("u", new ExprFieldWrapper(u0));
  w.addField("soln", new ExprFieldWrapper(soln));
  w.addField("u_x", new ExprFieldWrapper(du0_dx));
  w.addField("u_y", new ExprFieldWrapper(du0_dy));
  w.addField("u_z", new ExprFieldWrapper(du0_dz));
  w.write();

  Expr e0 = Integral(interior, soln, quad);
  Expr e1 = Integral(interior, du0_dx, quad);
  Expr e2 = Integral(interior, du0_dy, quad);
  Expr e3 = Integral(interior, du0_dz, quad);
  Expr I4 = Integral(interior, grad*F, quad);
  Expr I5 = Integral(bdry, F*n, quad);
  double err0 = fabs(evaluateIntegral(mesh, e0));
  double err1 = fabs(evaluateIntegral(mesh, e1));
  double err2 = fabs(evaluateIntegral(mesh, e2));
  double err3 = fabs(evaluateIntegral(mesh, e3));
  double val4 = fabs(evaluateIntegral(mesh, I4));
  double val5 = fabs(evaluateIntegral(mesh, I5));

  

  Out::os() << "error in Green's formula = " << err0 << std::endl;
  Out::os() << "error in du/dx = " << err1 << std::endl;
  Out::os() << "error in du/dy = " << err2 << std::endl;
  Out::os() << "error in du/dz = " << err3 << std::endl;

  Out::os() << "volume integral of div(F) = " << val4 << std::endl;
  Out::os() << "surface integral of F*n = " << val5 << std::endl;


    
  return SundanceGlobal::checkTest(err0 + err1 + err2 + err3, 1.0e-10);
}

/*
 * All done!
 */

#else



bool DiscFunc3D()
{
  std::cout << "dummy DiscFunc3D PASSED. Enable exodus to run the actual test" <<
    std::endl;
  return true;
}

#endif

