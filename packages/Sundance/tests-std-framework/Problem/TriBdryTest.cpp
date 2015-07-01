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
      Sundance::init(&argc, &argv);

      // Define vector type:
      VectorType<double> vecType = new EpetraVectorType();

      // Input Triangle mesh:
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource meshReader = new TriangleMeshReader("square.1", meshType );
      Mesh mesh = meshReader.getMesh();

      // Define mesh domain:
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter top = edges.labeledSubset(10);
      CellFilter right = edges.labeledSubset(20);
      CellFilter bottom = edges.labeledSubset(30);
      CellFilter left = edges.labeledSubset(40);

      CellFilter boundary = left + bottom + right + top; 

      Expr v = new TestFunction(new Lagrange(1));
      Expr u = new UnknownFunction(new Lagrange(1));
      Expr grad = gradient(2);

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr f = 2.0*(x*(1.0-x) + y*(1.0-y));

      QuadratureFamily quad = new GaussianQuadrature(2);

      Expr eqn = Integral(interior, (grad*v)*(grad*u) - f*v, quad);
      Expr bc = EssentialBC(boundary, v*u, quad);

      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

      LinearSolver<double> solver = LinearSolverBuilder::createSolver("amesos.xml");

      Expr soln = prob.solve(solver);

      Expr uEx = x*(1.0-x)*y*(1.0-y);

      FieldWriter w = new VTKWriter("soln");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(soln));
      w.write();

      double err = L2Norm(mesh, interior, soln-uEx, quad);
      
      Sundance::passFailTest(err, 1.0e-3);
    }
  catch(std::exception& e)
    {
      Sundance::handleException(e);
    }
  Sundance::finalize();
  return Sundance::testStatus();
}
