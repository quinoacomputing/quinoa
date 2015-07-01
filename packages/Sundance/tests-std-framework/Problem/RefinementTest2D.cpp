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
#include "SundanceUniformRefinementPair.hpp"




/* ----------------------------------------------------------------- */
//
// Demonstration and test of uniform refinement in 2D
//
/* ----------------------------------------------------------------- */






/* This function gets the DOF of a function at the specified node LID */
int getNodeDof(const RCP<DOFMapBase>& dofMap, int nodeLID)
{
  Array<int> dofs;
  dofMap->getDOFsForCell(0, nodeLID, 0, dofs);
  TEUCHOS_TEST_FOR_EXCEPT(dofs.size() != 1);
  return dofs[0];
}


int main(int argc, char** argv)
{
  
  try
  {
    Sundance::init(&argc, &argv);
    int np = MPIComm::world().getNProc();
    TEUCHOS_TEST_FOR_EXCEPT(np > 1); // works only in serial for now

    /* We will do our linear algebra using Epetra */
    VectorType<double> vecType = new EpetraVectorType();

    /* Create a mesh. This will be our coarsest mesh */
    MeshType meshType = new BasicSimplicialMeshType();
    int nx = 3;

    MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, 1,
      0.0, 1.0, nx, 1,
      meshType);
    Mesh coarse = mesher.getMesh();
    /* Label some random vertex. We'll check that the refinement preserves
     * vertex labels */
    coarse.setLabel(0, 3, 5);

    /* Refine the mesh. The URP object will store both meshes as well as
     * maps between cell indices at the two levels. */
    UniformRefinementPair urp(meshType, coarse);
    /* Run a consistency check on the refinement maps */
    int bad = urp.check();

    /* Get the fine mesh */
    const Mesh& fine = urp.fine();

    /* Create some expressions to be tested. The expression is linear
     * and so can be represented exactly with the P1 basis. Interpolation
     * to the fine mesh will also be exact. */
    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);
    Expr f = x + sqrt(2.0)*y;

    BasisFamily P1 = new Lagrange(1);
    
    /* Project a function onto the coarse space */
    DiscreteSpace coarseP1(coarse, P1, vecType);
    L2Projector projCoarseP1(coarseP1, f);
    Expr cP1 = projCoarseP1.project();
    Vector<double> cxP1 = DiscreteFunction::discFunc(cP1)->getVector();
    RCP<DOFMapBase> cDofMap = coarseP1.map();

    /* Project the same function onto the fine space */
    DiscreteSpace fineP1(fine, P1, vecType);
    L2Projector projFineP1(fineP1, f);
    Expr fP1 = projFineP1.project();
    Vector<double> fxP1 = DiscreteFunction::discFunc(fP1)->getVector();
    RCP<DOFMapBase> fDofMap = fineP1.map();

    /* Now we're going to interpolate the coarse function onto the fine
     * space. If the code is correct this will be equal to the function 
     * projected onto the fine mesh. */
    Vector<double> refinedX = fxP1.space().createMember();
    
    /* map the coarse space to the fine space */
    for (int f=0; f<fine.numCells(0); f++)
    {
      /* get the DOF for the node */
      int fDof = getNodeDof(fDofMap, f);
      /* If this fine vertex is on an edge of the coarse mesh, interpolate
       * from the facets of the edge. If the fine vertex is coincident
       * with a coarse vertex, use the value on that vertex */
      if (urp.newVertIsOnEdge()[f]) /* new vertex is on a coarse edge */
      {
        /* get the parent edge's LID in the coarse mesh */
        int cEdge = urp.newVertToOldLIDMap()[f];
        int ori; // dummy orientation, not needed here
        int cNode1 = coarse.facetLID(1, cEdge, 0, 0, ori);
        int cNode2 = coarse.facetLID(1, cEdge, 0, 1, ori);
        int cDof1 = getNodeDof(cDofMap, cNode1);
        int cDof2 = getNodeDof(cDofMap, cNode2);
        refinedX[fDof] = 0.5*(cxP1[cDof1]+cxP1[cDof2]);
      }
      else /* new vertex coincides with an old vertex */
      {
        int cNode = urp.newVertToOldLIDMap()[f];
        int cDof = getNodeDof(cDofMap, cNode);
        refinedX[fDof] = cxP1[cDof];
      }
    }
    
    /* Compare the R-mapped vector to the vector computed natively 
     * on the fine space */
    double err = (fxP1 - refinedX).norm2();
    cout << "error |fine - R*coarse| = " << err << endl;

    
    double finalTol = 0.5;
    Sundance::passFailTest(bad, finalTol);
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return Sundance::testStatus(); 
}
