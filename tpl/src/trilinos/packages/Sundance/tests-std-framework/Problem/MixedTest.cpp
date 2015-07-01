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

int main(int argc, char** argv)
{
  try
  {
    int nx = 1;
    int ny = 1;
    std::string solverFile = "amesos.xml";

    Sundance::init(&argc, &argv);
    int np = MPIComm::world().getNProc();

    nx = nx*np;
    ny = ny*np;

    VectorType<double> vecType = new EpetraVectorType();

    MeshType meshType = new BasicSimplicialMeshType();
      
    MeshSource mesher;

    int npx = -1;
    int npy = -1;
    PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
    TEUCHOS_TEST_FOR_EXCEPT(npx < 1);
    TEUCHOS_TEST_FOR_EXCEPT(npy < 1);
    TEUCHOS_TEST_FOR_EXCEPT(npx * npy != np);
    mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, npx, 
      0.0,  1.0, ny, npy, meshType);

    Mesh mesh = mesher.getMesh();

    CellFilter interior = new MaximalCellFilter();

    BasisFamily basis0 = new Lagrange(0);
    BasisFamily basis1 = new Lagrange(1);
    BasisFamily basis2 = new Lagrange(2);

    Expr ux = new UnknownFunction(basis0, "ux");
    Expr uy = new UnknownFunction(basis0, "uy");
    Expr vx=new TestFunction(basis0, "vx");
    Expr vy=new TestFunction(basis0, "vy");

    Expr sx = new UnknownFunction(basis1, "sx");
    Expr sy = new UnknownFunction(basis1, "sy");
    Expr tx=new TestFunction(basis1, "tx");
    Expr ty=new TestFunction(basis1, "ty");


    Expr p = new UnknownFunction(basis2, "p");
    Expr q = new TestFunction(basis2, "q");

    Expr dx = new Derivative(0);
    Expr dy = new Derivative(1);

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);

    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad = new GaussianQuadrature(4);

    /* Define the weak form */
    Expr exactP = 2.0*x*y + 3.0*x*x + 4.0*y*y;
    Expr exactS = List(2.0*y+6.0*x, 2.0*x+8.0*y);
    Expr exactU = List(6.0, 8.0);


    Expr eqn = Integral(interior,
      q*(p-exactP) + tx*(sx-dx*p) + ty*(sy-dy*p)
      + vx*(ux - dx*sx) + vy*(uy - dy*sy),
      quad);
    Expr bc;

    /* We can now set up the linear problem! */
    LinearProblem prob(mesh, eqn, bc, 
      List(q, vx, vy, tx, ty),
      List(p, ux, uy, sx, sy), vecType);

    Array<Expr> exact = tuple(exactP, exactU[0], exactU[1],
      exactS[0], exactS[1]);

    Out::os() << "row map" << std::endl;
    prob.rowMap(0)->print(Out::os());
    Out::os() << "col map" << std::endl;
    prob.colMap(0)->print(Out::os());


    ParameterXMLFileReader reader(solverFile);
    ParameterList solverParams = reader.getParameters();
    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver(solverParams);

    Expr soln = prob.solve(solver);

    FieldWriter w = new VTKWriter( "MixedTest2D" );
    w.addMesh( mesh );
    w.addField( "ux" , new ExprFieldWrapper( soln[1] ) );
    w.addField( "uy" , new ExprFieldWrapper( soln[2] ) );
    w.addField( "sx" , new ExprFieldWrapper( soln[3] ) );
    w.addField( "sy" , new ExprFieldWrapper( soln[4] ) );
    w.write();

    double totalErrSq = 0.0;
    for (int i=0; i<exact.size(); i++)
    {
      Expr err = exact[i] - soln[i];
      Expr errExpr = Integral(interior,err*err,quad);
      FunctionalEvaluator errInt(mesh, errExpr);
      double errorSq = errInt.evaluate();
      double err_i = std::sqrt(errorSq);
      Out::os() << "i=" << i << " error=" << err_i << std::endl;
      totalErrSq += errorSq;
    }

    double tol = 1.0e-12;
    Sundance::passFailTest(sqrt(totalErrSq), tol);

  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 

  return Sundance::testStatus();
}
