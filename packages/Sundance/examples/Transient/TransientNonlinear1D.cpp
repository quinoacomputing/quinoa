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

/* ------------------------------------------------------------------------ 
 *
 * This program solves a nonlinear initial-boundary-value problem in 
 * one spatial dimension using Galerkin finite elements in space and
 * Crank-Nicolson stepping in time. 
 *
 * ------------------------------------------------------------------------ */



Expr force(const double& epsilon, const Expr& x, const Expr& t)
{
  const double pi = 4.0*atan(1.0);

  Expr rtn = -8.0*epsilon*cos(2.0*pi*t)*
    pow(1.0 + pow(x,2.0)*epsilon*cos(2.0*pi*t),2.0)*
    (1.0 + 7.0*pow(x,2.0)*epsilon*cos(2.0*pi*t)) - 
    2.0*pi*pow(x,2.0)*epsilon*sin(2.0*pi*t);

  return rtn;
}


int main(int argc, char** argv)
{
  const double pi = 4.0*atan(1.0);

  try
  {
    Sundance::init(&argc, &argv);
    int np = MPIComm::world().getNProc();

    /* We will do our linear algebra using Epetra */
    VectorType<double> vecType = new EpetraVectorType();

    /* Create a mesh. It will be of type BasisSimplicialMesh, and will
     * be built using a PartitionedLineMesher. */
    MeshType meshType = new BasicSimplicialMeshType();
    int nx = 32;
    MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx*np, meshType);
    Mesh mesh = mesher.getMesh();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();
    CellFilter points = new DimensionalCellFilter(0);
    CellFilter leftPoint = points.coordSubset(0, 0.0);
    CellFilter rightPoint = points.coordSubset(0, 1.0);

      
    /* Create unknown and test functions, discretized using first-order
     * Lagrange interpolants */
    BasisFamily bas = new Lagrange(1);
    Expr u = new UnknownFunction(bas, "u");
    Expr v = new TestFunction(bas, "v");

    /* Create differential operator and coordinate function */
    Expr dx = new Derivative(0);
    Expr x = new CoordExpr(0);

    double epsilon = 0.25;
    /* The initial profile is u(x,0)=1 + epsilon x^2. 
     * Project this onto a discrete function. This discrete function,
     * uPrev, will represent
     * the initial conditions but will also be re-used as the starting value for
     * each timestep */
    Expr uStart = 1.0 + epsilon*x*x;
    DiscreteSpace discSpace(mesh, bas, vecType);
    L2Projector projector(discSpace, 1.0 + epsilon*x*x);
    Expr uPrev = projector.project();

    /* We need another discrete function for the current Newton approximation */
    Expr uNewt = copyDiscreteFunction(uPrev, "uNewt");

    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad = new GaussianQuadrature(4);

    int nSteps = nx;
    double dt = 1.0/((double) nSteps);
    /* Represent the time variable as a parameter expression, NOT as
     * a double variable. The reason is that we need to be able to update
     * the time value without rebuilding expressions. */
    Expr t = new Sundance::Parameter(0.0);
    Expr tPrev = new Sundance::Parameter(0.0);

    /* Define the weak form, semidiscretized in time */
    Expr eqn = Integral(interior, v*(u-uPrev) 
      + dt/2.0*(dx*v)*((dx*pow(u, 4.0))+(dx*pow(uPrev, 4.0)))
      - dt/2.0*v*(force(epsilon, x, t)+force(epsilon, x, tPrev)), quad); 
    /* Define the Dirichlet BC on the right side of the domain */
    Expr bc = EssentialBC(rightPoint, v*(u - 1.0 - epsilon*cos(2.0*pi*t)),quad);

    /* We can now set up the nonlinear problem! */
    NonlinearProblem prob(mesh, eqn, bc, v, u, uNewt, vecType); 

    NonlinearSolver<double> solver 
      = NonlinearSolverBuilder::createSolver("playa-newton-amesos.xml");

    /* Wrtie the initial conditions */
    FieldWriter w = new DSVWriter("transientNonlinear1D-0.dat"); 
    w.addMesh(mesh);
    w.addField("u", new ExprFieldWrapper(uPrev[0]));
    w.write();

    /* loop over timesteps */
    for (int i=0; i<nSteps; i++)
    {
      /* Set the times t_i and t_{i+1} */
      Out::root() << "timestep #" << i << endl;
      t.setParameterValue((i+1)*dt);
      tPrev.setParameterValue(i*dt);

      SolverState<double> state = prob.solve(solver);

      TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
        std::runtime_error,
        "Nonlinear solve failed to converge: message=" << state.finalMsg());
          
      updateDiscreteFunction(uNewt, uPrev);
      
      FieldWriter writer = new MatlabWriter("transientNonlinear1D-" 
        + Teuchos::toString(i+1) + ".dat");
      writer.addMesh(mesh);
      writer.addField("u", new ExprFieldWrapper(uPrev[0]));
      writer.write();
    }
    
    double err = L2Norm(mesh, interior, uPrev - uStart, quad);
    Out::root() << "dt=" << setw(16) << dt << " error = " << setw(16) << err << endl;
    
    double h = 1.0/(nx-1.0);
    double tol = 0.1*(dt*dt + h*h);
    Sundance::passFailTest(err, tol);
  }
  catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}
