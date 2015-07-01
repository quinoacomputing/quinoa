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

#if defined(HAVE_SUNDANCE_EXODUS) && defined(Trilinos_DATA_DIR)

int main( int argc , char **argv )
{
  try 
  {
    Sundance::init( &argc , &argv );
    
    VectorType<double> vecType = new EpetraVectorType();
    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new ExodusMeshReader("couette-2", meshType);
    Mesh mesh = mesher.getMesh();

    BasisFamily P2 = new Lagrange(2);
    BasisFamily P1 = new Lagrange(1);
    QuadratureFamily quad = new GaussianQuadrature(4);

    Expr ux = new UnknownFunction(P2);
    Expr uy = new UnknownFunction(P2);
    Expr vx = new TestFunction(P2);
    Expr vy = new TestFunction(P2);
    Expr u = List(ux, uy);
    Expr v = List(vx, vy);

    Expr p = new UnknownFunction(P1);
    Expr q = new TestFunction(P1);

    CellFilter interior = new MaximalCellFilter();
    CellFilter sides = new DimensionalCellFilter(1);
    CellFilter pts = new DimensionalCellFilter(0);
    CellFilter peg = pts.labeledSubset(1);
    CellFilter inner = sides.labeledSubset(2);
    CellFilter outer = sides.labeledSubset(1);

    Expr grad = gradient(2);
    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);
    Expr uInnerX = 0.0;
    Expr uInnerY = 0.0;
    Expr uOuterX = -y;
    Expr uOuterY = x;
    Expr uInner = List(uInnerX, uInnerY);
    Expr uOuter = List(uOuterX, uOuterY);

    DiscreteSpace uSpace(mesh, List(P2, P2), vecType);
    DiscreteSpace pSpace(mesh, P1, vecType);
    Expr uNewt = new DiscreteFunction(uSpace, 0.0, "u0");
    Expr pNewt = new DiscreteFunction(pSpace, 0.0, "p0");
    Expr upNewt = List(uNewt, pNewt);

    Expr nu = 1.0;
    double C1 = 25.0;
    double C2 = 25.0;

    /* Navier-Stokes equations */
    Expr eqn = Integral(interior, 
      v*((u*grad)*u)
      + nu * colonProduct( outerProduct(grad, u), outerProduct(grad, v) ) 
      - p*div(v) - q*div(u), quad)
      + NitscheStokesNoSlipBC(inner, quad, nu, v, q, u, p, uInner, C1, C2)
      + NitscheStokesNoSlipBC(outer, quad, nu, v, q, u, p, uOuter, C1, C2);

    Expr dummyBC;

    /* Equations for the PCD preconditioner */
    Expr MpEqn = Integral(interior, p*q, quad);
    Expr ApEqn = Integral(interior, (grad*p)*(grad*q), quad) 
      + Integral(peg, p*q, quad);
    Expr FpEqn = Integral(interior, q*((uNewt*grad)*p) 
      + (grad*p)*(grad*q), quad) 
      + Integral(peg, p*q, quad);

    /* Arrange the variables into blocks: [velocity, pressure] */
    Array<Block> unkBlocks = tuple(Block(u, vecType), Block(p, vecType));
    Array<Block> varBlocks = tuple(Block(v, vecType), Block(q, vecType));

    /* Form the PCD problems */
    LinearProblem MpProb(mesh, MpEqn, dummyBC, q, p, vecType);
    LinearProblem ApProb(mesh, ApEqn, dummyBC, q, p, vecType);
    LinearProblem FpProb(mesh, FpEqn, dummyBC, q, p, vecType);

    /* Read the parameters. This is a little complicated because the
     * PCD precond parameters contain solver parameter lists */
    ParameterXMLFileReader outerSolverReader("pcd-outer-belos.xml");
    ParameterList outerSolverParams = outerSolverReader.getParameters();
    ParameterList precParams = outerSolverParams.sublist("PCD Preconditioner");
    LinearSolver<double> outerSolver = LinearSolverBuilder::createSolver(outerSolverParams);

    /* Form the preconditioner, and set it as a user preconditioner */
    PreconditionerFactory<double> outerPrec 
      = new PCDPreconditionerFactory(precParams, MpProb, ApProb, FpProb);
    outerSolver.setUserPrec(outerPrec);

    /* set up the nonlinear problem */
    NonlinearProblem prob(mesh, eqn, dummyBC, varBlocks, unkBlocks,
      upNewt);
      

    int maxNewtIters = 10;
    double newtTol = 1.0e-6;
    bool newtonConverged = false;

    
    /* Allocate objects for the Jacobian, residual, and Newton step */
    LinearOperator<double> J = prob.allocateJacobian();
    Vector<double> resid = J.range().createMember();
    Vector<double> newtonStep = J.domain().createMember();

    for (int j=0; j<maxNewtIters; j++)
    {
      prob.setEvalPoint(upNewt);
      prob.computeJacobianAndFunction(J, resid);
      Out::root() << "|resid| = " << resid.norm2() << endl;
      if (resid.norm2() < newtTol) 
      {
        newtonConverged = true;
        break;
      }
      SolverState<double> solveState 
        = outerSolver.solve(J, -1.0*resid, newtonStep);

      TEUCHOS_TEST_FOR_EXCEPTION(solveState.finalState() != SolveConverged,
        std::runtime_error,
        "linear solve failed!");

      addVecToDiscreteFunction(upNewt, newtonStep);
      double newtStepNorm = newtonStep.norm2();
      Out::root() << "|newt step| = " << newtStepNorm << endl;
      if (newtStepNorm < newtTol) 
      {
        newtonConverged = true;
        break;
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(!newtonConverged, std::runtime_error,
      "Newton's method failed to converged after " << maxNewtIters
      << " iterations");

    Expr uSoln = uNewt;
    Expr pSoln = pNewt;
    

    /* Compare to the exact solution */
    Expr r = sqrt(x*x+y*y);
    Expr uMagExact = 4.0/3.0*r - 1.0/3.0/r;
    Expr uxExact = -uMagExact*y/r;
    Expr uyExact = uMagExact*x/r;
    Expr uExact = List(uxExact, uyExact);

    double errorNorm = L2Norm(mesh, interior, uSoln-uExact, quad);
  
    Out::root() << "error norm = " << errorNorm << std::endl;

    


    // project velocity onto P1
    DiscreteSpace P1Space( mesh , P1 , vecType );
    L2Projector proj_ux( P1Space , uSoln[0] );
    L2Projector proj_uy( P1Space , uSoln[1] );
    L2Projector proj_uxEx( P1Space , uxExact );
    L2Projector proj_uyEx( P1Space , uyExact );

    Expr ux_P1 = proj_ux.project( );
    Expr uy_P1 = proj_uy.project( );

    Expr uxEx_P1 = proj_uxEx.project( );
    Expr uyEx_P1 = proj_uyEx.project( );
    
    FieldWriter w = new VTKWriter("PCDNavierStokes");
    w.addMesh(mesh);
    w.addField( "ux" , new ExprFieldWrapper( ux_P1 ) );
    w.addField( "uy" , new ExprFieldWrapper( uy_P1 ) );
    w.addField( "p" , new ExprFieldWrapper( pSoln ) );
    w.addField( "uxEx" , new ExprFieldWrapper( uxEx_P1 ) );
    w.addField( "uyEx" , new ExprFieldWrapper( uyEx_P1 ) );
    w.write();

    Sundance::passFailTest(errorNorm, 0.001);
  }
  catch (std::exception &e) {
    Out::os() << e.what() << std::endl;
  }

  Sundance::finalize();
  return Sundance::testStatus(); 
}

    


#else


int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy PCDStokes PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif    
