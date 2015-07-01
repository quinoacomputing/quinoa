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


#include "SundanceExtrusionMeshTransformation.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;



Mesh ExtrusionMeshTransformation::apply(const Mesh& inputMesh) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(inputMesh.spatialDim() != 2, std::runtime_error,
                     "ExtrusionMeshTransformation::applyLocal() given mesh with "
                     "dimension " << inputMesh.spatialDim() << ". The "
                     "extrusion filter expects a 2D mesh as input");

  /* create the output 3D mesh, using the communicator from the input mesh */
  Mesh outputMesh = createMesh(3, inputMesh.comm());

  //  BasicSimplicialMesh* outputMeshP 
  //    = dynamic_cast<BasicSimplicialMesh*>(outputMesh.ptr().get());

	int nPts = inputMesh.numCells(0);
	
	/* create points in output mesh */
  int pointCount = 0;
  outputMesh.estimateNumVertices(nzLevels_ * nPts);

	for (int i=0; i<nPts; i++)
		{
			const Point& p = inputMesh.nodePosition(i);
			for (int j=0; j<=nzLevels_; j++, pointCount++)
				{
					double z = z0_ + (z1_-z0_)*j/((double) nzLevels_);
					outputMesh.addVertex(pointCount, Point(p[0], p[1], z), 0, 0);
				}
		}


  /* now do the elements */
	int nTri = inputMesh.numCells(2);

  /* each triangle in the input mesh will give rise to 
   * three tets per extrusion level in the output mesh */

  outputMesh.estimateNumElements(3*nTri*nzLevels_);

  Array<int> facetNodes(3);
  Array<int> facetOrientations(3);

  int tetCount = 0;

  TEUCHOS_TEST_FOR_EXCEPTION(inputMesh.cellType(2) != TriangleCell,
                     std::runtime_error,
                     "ExtrusionMeshTransformation::applyLocal() detected a "
                     "non-triangular mesh");
	for (int i=0; i<nTri; i++)
		{
      inputMesh.getFacetArray(2, i, 0, facetNodes, facetOrientations);

      SUNDANCE_OUT(this->verb()==3,
                   "triangle=" << i << " facet nodes are " << facetNodes);

      /* We sort the facets in order to get a consistent direction
       * for the edges at the different z-planes.  */
      std::sort(facetNodes.begin(), facetNodes.end());

      SUNDANCE_OUT(this->verb()==3,
                   "triangle=" << i << " sorted facet nodes are " 
                   << facetNodes);

			int a0 = facetNodes[0]*(nzLevels_+1);
			int b0 = facetNodes[1]*(nzLevels_+1);
			int c0 = facetNodes[2]*(nzLevels_+1);

			for (int j=0; j<nzLevels_; j++, tetCount+=3)
				{
					int a = a0 + j;
					int b = b0 + j;
					int c = c0 + j;
					int a1 = a+1;
					int b1 = b+1;
					int c1 = c+1;
          outputMesh.addElement(tetCount, tuple(a, a1, b1, c1), 0, 0);
          outputMesh.addElement(tetCount+1, tuple(a, b, b1, c1), 0, 0);
          outputMesh.addElement(tetCount+2, tuple(a, b, c, c1), 0, 0);
				}
		}

	return outputMesh;

}			

