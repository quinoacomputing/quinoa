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
#include "SundanceEvaluator.hpp"
#include "SundanceBubble.hpp"
#include "SundanceRivaraDriver.hpp"



int main(int argc, char** argv)
{
  
  try
    {
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh object to fill in from a file. It will be of type BasisSimplicialMesh. */
      MeshType meshType = new BasicSimplicialMeshType();

      /* read the mesh from the file disk*/
      MeshSource mesher = new ExodusMeshReader("disk", meshType);
      Mesh mesh = mesher.getMesh();

      Expr x = new CoordExpr(0);
      Expr err = exp(-2.0*x);

      DiscreteSpace discSpace(mesh, new Lagrange(0), vecType);
      L2Projector proj(discSpace, err);
      Expr errEst = proj.project();

      double goal = 1.0e-2;

      RefinementTransformation refiner(meshType, errEst, goal, 1.0e-4);
      mesh = refiner.apply(mesh);

      FieldWriter w = new VTKWriter("AMR");
      w.addMesh(mesh);
      w.write();

      double errNorm = 1.0;
      double tol = 1.0e-4;
      Sundance::passFailTest(errNorm, tol);
      
    }
  catch(std::exception& e)
    {
      std::cerr << e.what() << std::endl;
    }
  Sundance::finalize(); return Sundance::testStatus(); 
}


