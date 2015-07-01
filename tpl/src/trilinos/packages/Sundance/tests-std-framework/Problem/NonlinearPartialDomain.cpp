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
#include "SundanceCFMeshPair.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceInhomogeneousNodalDOFMap.hpp"


/** 
 * Tests logical operations on cell filters
 */


CELL_PREDICATE(ATest, {return x[0] <= 0.4;})

CELL_PREDICATE(BTest, {return x[0] >= 0.4 && x[0] <= 0.6;})

CELL_PREDICATE(CTest, {return x[0] >= 0.6;})
  
bool NonlinearPartialDomain()
{
  int np = MPIComm::world().getNProc();

  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Create a mesh. It will be of type BasisSimplicialMesh, and will
   * be built using a PartitionedLineMesher. */
  MeshType meshType = new BasicSimplicialMeshType();
  int nx = 200;
  int ny = 4;
  MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, np,
    0.0, 2.0, ny, 1,
    meshType);
  Mesh mesh = mesher.getMesh();

  Expr x = new CoordExpr(0);
  Expr dx = new Derivative(0);

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
  CellFilter bdry = new BoundaryCellFilter();
  CellFilter left = bdry.subset(new CoordinateValueCellPredicate(0,0.0));
  CellFilter right = bdry.subset(new CoordinateValueCellPredicate(0,1.0));

  CellFilter A = interior.subset(new ATest());
  CellFilter B = interior.subset(new BTest());
  CellFilter C = interior.subset(new CTest());
  CellFilter D = right;

  Expr u1 = new UnknownFunction(new Lagrange(1));
  Expr u2 = new UnknownFunction(new Lagrange(1));

  Expr v1 = new TestFunction(new Lagrange(1));
  Expr v2 = new TestFunction(new Lagrange(1));

  QuadratureFamily quad = new GaussianQuadrature(4);
  Expr eqn = Integral(interior, (dx*u1)*(dx*v1), quad) 
    + Integral(A, v2*(u2 - u1*u1), quad) 
    + Integral(B, v2*(u2 - 0.4*u1), quad) ;
  Expr bc = EssentialBC(left, v1*u1, quad) 
    + EssentialBC(right, v1*(u1-1.0), quad);

  /* Create a discrete space, and discretize the function 1.0 on it */
  Array<CellFilter> funcDomains = tuple(interior, A+B);
  BasisFamily L1 = new Lagrange(1);
  DiscreteSpace discSpace(mesh, BasisArray(tuple(L1, L1)), funcDomains, vecType);
  Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");



  /* We can now set up the nonlinear problem! */
  NonlinearProblem prob(mesh, eqn, bc, List(v1, v2), List(u1, u2), u0, vecType);


#ifdef HAVE_CONFIG_H
  ParameterXMLFileReader reader(searchForFile("SolverParameters/nox.xml"));
#else
  ParameterXMLFileReader reader("nox.xml");
#endif
  ParameterList noxParams = reader.getParameters();

  std::cerr << "solver params = " << noxParams << std::endl;

  NOXSolver solver(noxParams);

  prob.solve(solver);

  DiscreteSpace ds2(mesh, L1, A+B, vecType);
  L2Projector proj(ds2, 3.0*u0[1]);
  Expr uDisc = proj.project();

  Vector<double> vec = DiscreteFunction::discFunc(u0)->getVector();

  Expr err = Integral(interior, pow(u0[0] - x, 2.0), quad)
    + Integral(A, pow(u0[1] - u0[0]*u0[0], 2.0), quad)
    + Integral(B, pow(u0[1] - 0.4*u0[0], 2.0), quad);

  FunctionalEvaluator errInt(mesh, err);
  double errorSq = errInt.evaluate();
  std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

  Expr err2 = Integral(D, pow(dx*(u0[0] - x), 2.0), quad);
  FunctionalEvaluator err2Int(mesh, err);
  double error2Sq = err2Int.evaluate();
  std::cerr << "derivative error norm = " << sqrt(error2Sq) << std::endl << std::endl;

  /* Write the field in VTK format */
  FieldWriter w = new VTKWriter("PartialDomain2d");
  w.addMesh(mesh);
  w.addField("u1", new ExprFieldWrapper(u0[0]));
  w.addField("u2", new ExprFieldWrapper(u0[1]));
  w.write();

  return SundanceGlobal::checkTest(sqrt(errorSq), 1.0e-4);
}

