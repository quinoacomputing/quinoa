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


#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceExodusMeshReader.hpp"
#include "SundanceExodusWriter.hpp"
#include "SundanceTriangleWriter.hpp"
#include "SundanceVTKWriter.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;



static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}


class ProcField : public FieldBase
{
public: 
  /** */
  ProcField(int dim)
    : dim_(dim), pid_(MPIComm::world().getRank())
    {}

  /** */
  virtual ~ProcField() {;}

  /** */
  virtual double getData(int cellDim, int cellID, int elem) const 
    {
      TEUCHOS_TEST_FOR_EXCEPT(!isDefined(cellDim, cellID, elem));
      return pid_;
    }
  
  /** */
  virtual bool isDefined(int cellDim, int cellID, int elem) const 
    {
      return cellDim == dim_;
    }
  
  /** */
  virtual bool isPointData() const {return dim_ == 0 ;}

  /* */
  GET_RCP(FieldBase);

private:
  int dim_;
  int pid_;
};

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusMeshReader("plate3D-0", meshType);



      Mesh mesh = mesher.getMesh();

      FieldWriter w = new VTKWriter("wheel");
      w.addMesh(mesh);
      w.addField("cell proc", new ProcField(mesh.spatialDim()));
      w.addField("point proc", new ProcField(0));
      w.write();

      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
      std::cerr << "Detected exception: " << e.what() << std::endl;
		}
}

