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

#ifndef SUNDANCE_GEOMUTILS_H
#define SUNDANCE_GEOMUTILS_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include <list>

namespace Sundance
{
  /** \relates Mesh 
   * Return the LID of the maximal cell that contains the point x.
   * If the point lays on an edge between two or more
   * cells, this will return the first of the those cells encountered
   * on a breadth-first search starting at the initial guess.
   */
  int findEnclosingCell(const Mesh& mesh, 
			int cellDim,
			int initialGuessLID, 
			const double* x);

  /** \relates Mesh
   * Test whether a point is enclosed in a cell. The cell is assumed
   * to be simplicial.
   */
  bool cellContainsPoint(const Mesh& mesh, 
			 int cellDim,
			 int cellLID, const double* x,
			 Array<int>& facetLID);

  /** \relates Mesh
   * Tests orientation of a point relative to a line.
   */
  double orient2D(const double* a, const double* b, const double* x);

  /** \relates Mesh 
   * Get the list of maximal neighbors of a cell 
   */
  void maximalNeighbors(const Mesh& mesh, int cellDim,
			int cellLID, const Array<int>& facetLID,
			std::list<int>& rtn);

  /** \relates Mesh 
   * Pullback a point to local coordinates within a cell
   */
  Point pullback(const Mesh& mesh, int cellDim, int cellLID, const double* x);

  /** */
  void printCell(const Mesh& mesh, int cellLID);

  /** */
  double volume(const Mesh& mesh, int cellDim, int cellLID);
  
  /** \relates Mesh */
  int lookupEdgeLIDFromVerts(const Mesh& mesh, int v1, int v2);
}


#endif


