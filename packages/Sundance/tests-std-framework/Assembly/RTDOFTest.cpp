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
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceRaviartThomas.hpp"
#include "SundanceFunctionSupportResolver.hpp"


int main(int argc, char** argv)
{
  int stat = 0;
  try
  {
    Sundance::init(&argc, &argv);

    /* */
    int npx = MPIComm::world().getNProc();

    /* Create a mesh. It will be of type BasisSimplicialMesh, and will
     * be built using a PartitionedLineMesher. */
    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedRectangleMesher(0.0, 3.0, 2, npx,
      0.0, 1.0, 2, 1,
      meshType);
    Mesh mesh = mesher.getMesh();


    FieldWriter w = new VerboseFieldWriter();
    w.addMesh(mesh);
    w.write();

    FieldWriter w2 = new VTKWriter("rtMesh");
    w2.addMesh(mesh);
    w2.write();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();

    /* */
    Expr v = new TestFunction(new RaviartThomas(2));
    /* */
    Expr u = new UnknownFunction(new RaviartThomas(2));

    /* */
    Expr q = new TestFunction(new Lagrange(1));
    /* */
    Expr p = new UnknownFunction(new Lagrange(1));

    Out::os() << "v = " << describeFunction(v) << std::endl;
    Out::os() << "u = " << describeFunction(u) << std::endl;

    Out::os() << "q = " << describeFunction(q) << std::endl;
    Out::os() << "p = " << describeFunction(p) << std::endl;



    QuadratureFamily quad = new GaussianQuadrature(2);
    Expr eqn = Integral(interior, v*u + p*q, quad);
    Expr dum;

    RCP<FunctionSupportResolver> fsr 
      = rcp(new FunctionSupportResolver(eqn, dum, tuple(List(v,q).flatten()), tuple(List(u,p).flatten()),
          dum, dum, tuple(dum), false));
          
    
    int verb = 0;
    DOFMapBuilder builder(mesh, fsr, false, verb);

    for (int br = 0; br<builder.rowMap().size(); br++)
    {
      RCP<DOFMapBase> rm = builder.rowMap()[br];
      rm->print(Out::os());
    }

    
  }
	catch(std::exception& e)
  {
    stat = -1;
    std::cerr << "RT dof test FAILED" << std::endl;
    std::cerr << e.what() << std::endl;
  }

  return stat;
  
}
