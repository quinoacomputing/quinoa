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
    const double pi = 4.0*atan(1.0);
    double lambda = 1.25*pi*pi;

    int nx = 32;
    int nt = 10;
    double tFinal = 1.0/lambda;

    Sundance::setOption("nx", nx, "Number of elements");
    Sundance::setOption("nt", nt, "Number of timesteps");
    Sundance::setOption("tFinal", tFinal, "Final time");
    
    Sundance::init(&argc, &argv);

    /* Creation of vector type */
    VectorType<double> vecType = new EpetraVectorType();

    /* Set up mesh */
    MeshType meshType = new BasicSimplicialMeshType();
      
    MeshSource meshSrc = new PartitionedRectangleMesher(
      0.0, 1.0, nx,
      0.0, 1.0, nx,
      meshType);
    Mesh mesh = meshSrc.getMesh();

    /* 
     * Specification of cell filters
     */
    CellFilter interior = new MaximalCellFilter();
    CellFilter edges = new DimensionalCellFilter(1);
    CellFilter west = edges.coordSubset(0, 0.0);
    CellFilter east = edges.coordSubset(0, 1.0);
    CellFilter south = edges.coordSubset(1, 0.0);
    CellFilter north = edges.coordSubset(1, 1.0);

    /* set up test and unknown functions */
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "u");
    Expr v = new TestFunction(basis, "v");

    /* set up differential operators */
    Expr grad = gradient(2);

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);

    Expr t = new Sundance::Parameter(0.0);
    Expr tPrev = new Sundance::Parameter(0.0);


    DiscreteSpace discSpace(mesh, basis, vecType);
    Expr uExact = cos(0.5*pi*y)*sin(pi*x)*exp(-lambda*t);
    L2Projector proj(discSpace, uExact);
    Expr uPrev = proj.project();


    /* 
     * We need a quadrature rule for doing the integrations 
     */
    QuadratureFamily quad = new GaussianQuadrature(2);

    double deltaT = tFinal/nt;

    Expr gWest = -pi*exp(-lambda*t)*cos(0.5*pi*y);
    Expr gWestPrev = -pi*exp(-lambda*tPrev)*cos(0.5*pi*y);
    
    /* Create the weak form */
    Expr eqn = Integral(interior, v*(u-uPrev)/deltaT
      + 0.5*(grad*v)*(grad*u + grad*uPrev), quad)
      + Integral(west, -0.5*v*(gWest+gWestPrev), quad);

    Expr bc = EssentialBC(east + north, v*u, quad);

    
    LinearProblem prob(mesh, eqn, bc, v, u, vecType);

    
    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver("amesos.xml");

    FieldWriter w0 = new VTKWriter("TransientHeat2D-0");
    w0.addMesh(mesh);
    w0.addField("T", new ExprFieldWrapper(uPrev[0]));
    w0.write();

    for (int i=0; i<nt; i++)
    {
      t.setParameterValue((i+1)*deltaT);
      tPrev.setParameterValue(i*deltaT);
      Out::root() << "t=" << (i+1)*deltaT << endl;
      Expr uNext = prob.solve(solver);
      
      std::ostringstream oss;
      oss << "TransientHeat2D-" << i+1;
      FieldWriter w = new VTKWriter(oss.str());
      w.addMesh(mesh);
      w.addField("T", new ExprFieldWrapper(uNext[0]));
      w.write();

      updateDiscreteFunction(uNext, uPrev);
    }


    
    double err = L2Norm(mesh, interior, uExact-uPrev, quad);
    Out::root() << "error norm=" << err << endl;

    double h = 1.0/(nx-1.0);
    double tol = 0.1*(pow(h,2.0) + pow(lambda*deltaT, 2.0));
    Out::root() << "tol=" << tol << endl;
    
    
    Sundance::passFailTest(err, tol);
  }
	catch(std::exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return Sundance::testStatus();
}

