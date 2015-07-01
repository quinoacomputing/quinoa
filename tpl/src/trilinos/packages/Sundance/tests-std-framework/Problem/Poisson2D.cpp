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


/* 
 * <Ignore> 
 * 
 *</Ignore>
 */

/* <Ignore> */
#include "Sundance.hpp"
using Sundance::List;

#if defined(HAVE_SUNDANCE_EXODUS)

/* </Ignore> */

/*
 * \documentclass[10pt]{article}
 * \begin{document}
 * 
 */


/* 
 * Solves the Poisson equation in 2D
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;}) 
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;}) 
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;}) 



int main(int argc, char** argv)
{
  try
  {
    int nx = 2;
    int ny = 2;
    std::string meshFile="builtin";
    std::string solverFile = "aztec-ml.xml";
    Sundance::setOption("meshFile", meshFile, "mesh file");
    Sundance::setOption("nx", nx, "number of elements in x");
    Sundance::setOption("ny", ny, "number of elements in y");
    Sundance::setOption("solver", solverFile, "name of XML file for solver");

    Sundance::init(&argc, &argv);
    int np = MPIComm::world().getNProc();

    nx = nx*np;
    ny = ny*np;

    /* 
     * <Header level="subsubsection" name="vector_type">
     * Creation of vector type
     * </Header>
     * 
     * Next we create a \verb+VectorType+ object. A \verb+VectorType+
     * is an abstract factory which produces \verb+VectorSpace+ objects.
     * A \verb+VectorSpace+ is itself a factory which can produce 
     * \verb+Vector+ objects which are, logically enough, vectors.
     * 
     * Why this heirarchy of factorys? Vectors may need to be created
     * inside code far below the user level, and the \verb+VectorSpace+ 
     * encapsulates the information needed to build a vector of a 
     * given type, size, and distribution over processors. The 
     * \verb+VectorType+ is needed because we might need to select a 
     * particular {\it software representation} for a vector, {\it e.g.}, 
     * Trilinos, Petsc, or some other library. 
     *
     * By using an \verb+EpetraVectorType+, we will be creating
     * \verb+EpetraVectorSpace+ vector spaces, which in turn
     * create \verb+EpetraVector+ vectors.
     */
    VectorType<double> vecType = new EpetraVectorType();

    /* 
     * <Header level="subsubsection" name="mesh">
     * Creation of mesh
     * </Header>
     * 
     * The creation of a mesh object involves several intermediate
     * objects: a \verb+MeshType+ and a \verb+MeshSource+. The 
     * \verb+MeshType+ specifies what mesh data structure will be used,
     * and the \verb+MeshSource+ builds that data structure. How it is
     * built will depend on the \verb+MeshSource+ subtype used. An
     * \verb+ExodusMeshReader+, for instance, reads from an Exodus II file,
     * while a \verb+PartitionedRectangleMesher+ builds a uniform 
     * triangulation of a rectangle.
     * 
     * Once the \verb+MeshType+ and \verb+MeshSource+ objects have been
     * defined, the mesh is created (read, built, whatever) by a call
     * to the \verb+getMesh()+ method of \verb+MeshSource+.
     */
    MeshType meshType = new BasicSimplicialMeshType();
      
    MeshSource mesher;
    if (meshFile != "builtin")
    {
      mesher = new ExodusMeshReader("../../../examples-tutorial/meshes/"+meshFile, meshType);
    }
    else
    {
      int npx = -1;
      int npy = -1;
      PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
      TEUCHOS_TEST_FOR_EXCEPT(npx < 1);
      TEUCHOS_TEST_FOR_EXCEPT(npy < 1);
      TEUCHOS_TEST_FOR_EXCEPT(npx * npy != np);
      mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, npx, 
        0.0,  1.0, ny, npy, meshType);
    }
    Mesh mesh = mesher.getMesh();


    bool meshOK = mesh.checkConsistency(meshFile+"-check");
    if (meshOK) 
    {
      cout << "mesh is OK" << std::endl;
    }
    else
    {
      cout << "mesh is INCONSISTENT" << std::endl;
    }
    mesh.dump(meshFile+"-dump");

    WatchFlag watchMe("watch eqn");
    watchMe.setParam("integration setup", 6);
    watchMe.setParam("integration", 6);
    watchMe.setParam("fill", 6);
    watchMe.setParam("evaluation", 0);
    watchMe.deactivate();

    WatchFlag watchBC("watch BCs");
    watchBC.setParam("integration setup", 6);
    watchBC.setParam("integration", 6);
    watchBC.setParam("fill", 6);
    watchBC.setParam("evaluation", 0);
//    watchBC.deactivate();
    


    /* 
     * <Header level="subsubsection" name="cell_filter">
     * Specification of geometric regions
     * </Header>
     * 
     * We'll need to specify subsets of the mesh on which equations
     * or boundary conditions are defined. In many FEA codes this is
     * done by explicit definition of element blocks, node sets, and
     * side sets. Rather than working with sets explicitly at the user
     * level, we instead work with filtering rules that produce 
     * sets of cells. These rules are represented by \verb+CellFilter+
     * objects. 
     *
     * \verb+MaximalCellFilter+
     * selects all cells having dimension equal to the spatial dimension
     * of the mesh. 
     */
    CellFilter interior = new MaximalCellFilter();

    
    /* 
     * \verb+DimensionalCellFilter+
     * selects all cells of a specified dimension.
     */
    CellFilter edges = new DimensionalCellFilter(1);

    CellFilter left;
    CellFilter right;
    CellFilter top;
    CellFilter bottom;

    if (meshFile != "builtin")
    {
      left = edges.labeledSubset(1);
      right = edges.labeledSubset(2);
      top = edges.labeledSubset(3);
      bottom = edges.labeledSubset(4);
    }
    else
    {
      left = edges.subset(new LeftPointTest());
      right = edges.subset(new RightPointTest());
      top = edges.subset(new TopPointTest());
      bottom = edges.subset(new BottomPointTest());
    }
      
    /* Create unknown and test functions, discretized using second-order
     * Lagrange interpolants */
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "u");
    Expr v = new TestFunction(basis, "v");

    /* Create differential operator and coordinate functions */
    Expr dx = new Derivative(0);
    Expr dy = new Derivative(1);
    Expr grad = List(dx, dy);
    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);

    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad2 = new GaussianQuadrature(2);
    QuadratureFamily quad4 = new GaussianQuadrature(4);

    /* Define the weak form */
    //Expr eqn = Integral(interior, (grad*v)*(grad*u) + v, quad);
    Expr one = new Sundance::Parameter(1.0);
    Expr exactSoln = 1.0 + sqrt(3.0)*y+sqrt(2.0)*x;
    Expr eqn = Integral(interior, (grad*u)*(grad*v), quad2, watchMe);
    /* Define the Dirichlet BC */
    Expr h = new CellDiameterExpr();
    Expr bc = EssentialBC(bottom+top+left+right, v*(u-exactSoln), quad4, watchBC);

    /* We can now set up the linear problem! */
    LinearProblem prob(mesh, eqn, bc, v, u, vecType);

    Out::os() << "row map" << endl;
    prob.rowMap(0)->print(Out::os());
    Out::os() << endl << "col map" << endl;
    prob.colMap(0)->print(Out::os());


    ParameterXMLFileReader reader(solverFile);
    ParameterList solverParams = reader.getParameters();
    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver(solverParams);

    Expr soln = prob.solve(solver);

    DiscreteSpace discSpace2(mesh, new Lagrange(2), vecType);
    DiscreteSpace discSpace0(mesh, new Lagrange(0), vecType);
    L2Projector proj1(discSpace2, soln-exactSoln);
    double pid = MPIComm::world().getRank();
    L2Projector proj2(discSpace0, pid);
    Expr errorDisc = proj1.project();
    Expr pidDisc = proj2.project();




    /* Write the field in VTK format */
    FieldWriter w = new VTKWriter("Poisson2d");
    w.addMesh(mesh);
    w.addField("soln", new ExprFieldWrapper(soln[0]));
    w.addField("error", new ExprFieldWrapper(errorDisc));
    w.addField("rank", new ExprFieldWrapper(pidDisc));
    w.write();

    FieldWriter w2 = new VerboseFieldWriter("mesh");
    w2.addMesh(mesh);
    w2.write();

    Expr err = exactSoln - soln;
    Expr errExpr = Integral(interior, 
      err*err,
      quad4);

    Expr derivErr = dx*(exactSoln-soln);
    Expr derivErrExpr = Integral(interior, 
      derivErr*derivErr, 
      quad2);

      

    Expr fluxErrExpr = Integral(top, 
      pow(dy*(soln-exactSoln), 2),
      new GaussianQuadrature(2));

      
//    watchBC.deactivate();
    Expr exactFluxExpr = Integral(top, 
      dy*exactSoln,
      new GaussianQuadrature(2), watchBC);

    Expr numFluxExpr = Integral(top, 
      dy*soln,
      new GaussianQuadrature(2));


    FunctionalEvaluator errInt(mesh, errExpr);
    FunctionalEvaluator derivErrInt(mesh, derivErrExpr);

    double errorSq = errInt.evaluate();
    cout << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

    double derivErrorSq = derivErrInt.evaluate();
    cout << "deriv error norm = " << sqrt(derivErrorSq) << std::endl << std::endl;


    double fluxErrorSq = evaluateIntegral(mesh, fluxErrExpr);
    cout << "flux error norm = " << sqrt(fluxErrorSq) << std::endl << std::endl;


    cout << "exact flux = " << evaluateIntegral(mesh, exactFluxExpr) << std::endl;
    cout << "numerical flux = " << evaluateIntegral(mesh, numFluxExpr) << std::endl;

    Sundance::passFailTest(sqrt(errorSq + derivErrorSq + fluxErrorSq), 1.0e-9);

  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 

  return Sundance::testStatus();
}


/* <Ignore> */
#else


int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy Poisson2D PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif

/* </Ignore> */

/* \end{document} */
