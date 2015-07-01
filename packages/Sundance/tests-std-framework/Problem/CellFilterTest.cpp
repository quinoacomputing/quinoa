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


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;});

CELL_PREDICATE(ATest, {return x[0] <= 0.4;});

CELL_PREDICATE(BTest, {return x[0] >= 0.4 && x[0] <= 0.6;});

CELL_PREDICATE(CTest, {return x[0] >= 0.6;});
  
int main(int argc, void** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 10;
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();

      Expr x = new CoordExpr(0);

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());
      CellFilter A = interior.subset(new ATest());
      CellFilter B = interior.subset(new BTest());
      CellFilter C = interior.subset(new CTest());

      Expr u1 = new UnknownFunction(new Lagrange(1));
      Expr u2 = new UnknownFunction(new Lagrange(1));

      Expr v1 = new TestFunction(new Lagrange(1));
      Expr v2 = new TestFunction(new Lagrange(1));

      QuadratureFamily quad = new GaussianQuadrature(2);
      Expr eqn = Integral(interior, v1*(u1 - 2.0), quad) 
        + Integral(A, v2*(u2 - x*u1), quad) 
        + Integral(B, v2*(u2 - 0.4*u1), quad) ;
      Expr bc;

      LinearProblem prob(mesh, eqn, bc, List(v1, v2), List(u1, u2), vecType);


      cout << "-------- ROW MAP --------------------------" << std::endl;
      prob.rowMap(0)->print(cout);

      cout << "-------- COL MAP --------------------------" << std::endl;
      prob.colMap(0)->print(cout);

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/aztec.xml"));
#else
      ParameterXMLFileReader reader("aztec.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      cout << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      cout << "matrix = " << std::endl << prob.getOperator() << std::endl;
      cout << "RHS = " << std::endl << prob.getRHS() << std::endl;

      Expr soln = prob.solve(solver);

      Vector<double> vec = DiscreteFunction::discFunc(soln)->getVector();

      cout << "solution = " << vec << std::endl;
      
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
