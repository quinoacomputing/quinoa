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

/* 
 * test for large numbers of cell filters
 */

int main(int argc, char** argv)
{
  try
  {
    int nx = 16;
    Sundance::setOption("nx", nx, "Number of elements");

    Sundance::init(&argc, &argv);

    VectorType<double> vecType = new EpetraVectorType();

    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
    Mesh mesh = mesher.getMesh();

    CellFilter interior = new MaximalCellFilter();
    CellFilter points = new DimensionalCellFilter(0);
    
    Out::os() << "making vert filters" << endl;
    Array<CellFilter> X(nx+1);
    for (int i=0; i<=nx; i++)
    {
      X[i] = points.coordSubset(0, i/((double) nx));
    }
    Out::os() << "combining vert filters" << endl;
    CellFilter nodes = X[0];
    for (int i=0; i<nx; i++) 
    {
      nodes = nodes + X[i+1];
    }
    
    
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "w");
    Expr v = new TestFunction(basis, "v");

    Expr x = new CoordExpr(0);

    QuadratureFamily quad2 = new GaussianQuadrature(2);

    Expr eqn 
      = Integral(interior, v*(u-x), quad2) + Integral(nodes, v*(u-x), quad2);

    Expr bc;

    Out::os() << "forming prob" << endl;

    LinearProblem prob(mesh, eqn, bc, v, u, vecType);

    LinearSolver<double> linSolver 
      = LinearSolverBuilder::createSolver("amesos.xml");

    Out::os() << "solving prob" << endl;
    Expr soln = prob.solve(linSolver);
    Out::os() << "done solving prob" << endl;
  }
	catch(std::exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return 0;
}

