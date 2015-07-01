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


#include "SundancePartitionedLineMesher.hpp"
#include "SundanceOut.hpp"
#include "SundanceCollectiveExceptionCheck.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;


PartitionedLineMesher::PartitionedLineMesher(const ParameterList& params)
  : MeshSourceBase(params),
    ax_(params.get<double>("ax")),
    bx_(params.get<double>("bx")),
    nx_(params.get<int>("nx"))
{;}

Mesh PartitionedLineMesher::fillMesh() const
{
  Mesh mesh;

  try
    {
      SUNDANCE_OUT(this->verb() > 0,
                   "PartitionedLineMesher::fillLocalMesh() is meshing "
                   "interval [" << ax_ << ", " << bx_ << "]");

      mesh = createMesh(1);


      /* compute number of points per proc */

      int np = nProc();
      int nppx = nx_/np;

      SUNDANCE_OUT(this->verb() > 0,
                   "PartitionedLineMesher::fillLocalMesh() has " << nppx
                   << " points per proc");

      int px = myRank();

      int lowestVisiblePtX = px*nppx-1;
      if (lowestVisiblePtX < 0) lowestVisiblePtX = 0;
	
      int highestVisiblePtX = lowestVisiblePtX + nppx + 1;
      if (highestVisiblePtX > nx_) highestVisiblePtX = nx_;

      SUNDANCE_OUT(this->verb() > 0,
                   "index range is [" << lowestVisiblePtX << ", " << 
                   highestVisiblePtX << "]");

      Array<int> pts(highestVisiblePtX-lowestVisiblePtX+1); 
      int globalIndex = 0;

      /* add the visible points into the mesh */
      for (int i=0; i<=nx_; i++, globalIndex++)
        {
          if (i < lowestVisiblePtX || i > highestVisiblePtX) continue;
          int pointOwner = i/nppx;
          if (i==nx_) pointOwner--;
          Point x( ax_ + ((double) i)*(bx_-ax_)/((double) nx_)); 

          SUNDANCE_OUT(this->verb() > 1, "adding point GID=" 
                       << globalIndex << " x=" << x << " owner=" << pointOwner); 
          int lid = mesh.addVertex(globalIndex, x, pointOwner, 0);
          pts[i-lowestVisiblePtX] = globalIndex;
          SUNDANCE_OUT(this->verb() >  3,
                       "point " << x << " registered with LID=" << lid);
        }

      /* add the visible cells to the mesh */
      globalIndex = 0 ;

      for (int i=0; i<nx_; i++, globalIndex++)
        {
          if (i < lowestVisiblePtX || i >= highestVisiblePtX) continue;
          int a = pts[i-lowestVisiblePtX];
          int b = pts[i-lowestVisiblePtX+1];
          int cellOwner = i/nppx;
          SUNDANCE_OUT(this->verb() > 1, "adding elem GID=" 
                       << globalIndex << " nodes=" << tuple(a,b) 
                       << " owner=" << cellOwner); 

          int lid = mesh.addElement(globalIndex, tuple(a,b), cellOwner, 0);
          SUNDANCE_OUT(this->verb() >  3,
                       "elem " << tuple(a,b) << " registered with LID=" << lid);
        }

      mesh.freezeTopology();

      if (px==0) mesh.setLabel(0, 0, 1); 
      if (px==np-1) mesh.setLabel(0, mesh.mapGIDToLID(0, nx_), 2);
    
  

    }
  catch(std::exception& e0)
    {
      reportFailure(comm());
      SUNDANCE_TRACE_MSG(e0, "while meshing a line");
    }
  TEUCHOS_TEST_FOR_EXCEPTION(checkForFailures(comm()), std::runtime_error, 
                     "off-proc error detected on proc=" << myRank()
                     << " while meshing line");
  return mesh;
  
}
