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

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})




int main( int argc , char **argv )
{
  try {
    
    Sundance::init( &argc , &argv );


#ifdef HAVE_CONFIG_H
    ParameterXMLFileReader reader(searchForFile("SolverParameters/anasazi-ml.xml"));
#else
    ParameterXMLFileReader reader("anasazi-ml.xml");
#endif

    ParameterList params = reader.getParameters();


    int np = MPIComm::world().getNProc();
    int npx = -1;
    int npy = -1;
    PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
    TEUCHOS_TEST_FOR_EXCEPT(npx < 1);
    TEUCHOS_TEST_FOR_EXCEPT(npy < 1);
    TEUCHOS_TEST_FOR_EXCEPT(npx * npy != np);

    Out::root() << "npx=" << npx << ", npy=" << npy << std::endl;

    VectorType<double> vecType = new EpetraVectorType();

    const int nx = 16;
    const int ny = 16;
    const double C = 2.0;
    bool lumpedMass = true;
    ParameterList solverParams = params.sublist("Eigensolver");

    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedRectangleMesher( 0.0 , 1.0 , nx , npx ,
      0.0 , 1.0 , ny , npy,
      meshType );
    Mesh mesh = mesher.getMesh();
    int dim = mesh.spatialDim();

    CellFilter interior = new MaximalCellFilter();
    CellFilter boundary = new BoundaryCellFilter();
    CellFilter left = boundary.subset( new LeftPointTest() );
    CellFilter right = boundary.subset( new RightPointTest() );
    CellFilter top = boundary.subset( new TopPointTest() );
    CellFilter bottom = boundary.subset( new BottomPointTest() );
    

    BasisFamily L = new Lagrange(1);
    Expr u = new UnknownFunction(L , "u" );
    Expr v = new TestFunction(L , "v" );
    QuadratureFamily quad = new GaussianQuadrature( 2 );
    
    Expr h = new CellDiameterExpr();
    Expr alpha = 4.0* C /h; 
    Expr n = CellNormalExpr(dim, "n");
    Expr grad = gradient(dim);

    
    /* weak form for the stiffness part of the Helmholtz equation */
    Expr eqn = Integral(interior, (grad*v)*(grad*u), quad)
      + Integral(left+right+top+bottom, 
        alpha*u*v - v*(n*grad)*u - u*(n*grad)*v, quad);

    LinearEigenproblem prob(mesh, eqn, v, u, vecType, lumpedMass);
    
    Eigensolver<double> solver = new AnasaziEigensolver<double>(solverParams);
    Eigensolution soln = prob.solve(solver);

    FieldWriter w = new VTKWriter( "Eigen2D" );
    w.addMesh( mesh );
    for (int i=0; i< soln.numEigenfunctions(); i++)
    {
      Expr ev = soln.eigenfunction(i);
      const std::complex<double>& ew = soln.eigenvalue(i);
      Out::root() << "ew=(" << ew  << ")" << std::endl;
      w.addField("u[" + Teuchos::toString(i) + "]", 
        new ExprFieldWrapper(ev));
    }
    w.write();

    const double pi = 4.0*atan(1.0);
    double exactEW = 2.0*pi*pi;

    double error = std::abs(exactEW - soln.eigenvalue(0))/std::abs(exactEW + soln.eigenvalue(0));
    
    double tol = 0.001*::pow(32/((double) nx), 2.0);

    Sundance::passFailTest(error, tol);
  }
  catch (std::exception &e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return Sundance::testStatus();
}

