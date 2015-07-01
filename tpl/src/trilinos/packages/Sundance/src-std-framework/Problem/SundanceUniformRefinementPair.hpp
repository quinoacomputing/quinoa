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

#ifndef SUNDANCE_UNIFORM_REFINEMENT_PAIR_H
#define SUNDANCE_UNIFORM_REFINEMENT_PAIR_H

#include "SundanceDefs.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceArrayOfTuples.hpp"

namespace Sundance
{
/**
 *
 */
class UniformRefinementPair
{
public:
  /** */
  UniformRefinementPair();
  /** */
  UniformRefinementPair(const MeshType& meshType,
    const Mesh& coarse);

  /** */
  const Mesh& fine() const {return fine_;}

  /** */
  const Mesh& coarse() const {return coarse_;}

  /** */
  const Array<int>& oldToNewVertMap() const {return oldToNewVertMap_;}

  /** */
  const Array<int>& newVertToOldLIDMap() const {return newVertToOldLIDMap_;}

  /** */
  const Array<int>& newVertIsOnEdge() const {return newVertIsOnEdge_;}

  /** */
  const Array<int>& oldEdgeToNewVertMap() const {return oldEdgeToNewVertMap_;}

  /** */
  const ArrayOfTuples<int>& oldToNewElemMap() const 
    {return oldToNewElemMap_;}

  /** */
  const Array<int>& newToOldElemMap() const 
    {return newToOldElemMap_;}

  /** */
  const Array<Array<int> >& oldEdgeChildren() const 
    {return oldEdgeChildren_;}

  /** */
  const Array<Array<int> >& oldEdgeParallels() const 
    {return oldEdgeParallels_;}
  /** */
  const Array<int>& newEdgeParents() const 
    {return newEdgeParents_;}

  /** */
  const Array<int>& newEdgeParallels() const 
    {return newEdgeParallels_;}


  /** */
  const ArrayOfTuples<int>& interiorEdgesOfCoarseElems() const 
    {return interiorEdges_;}


  /** Run a consistency check on the pair of meshes. Returns the number
   * of errors detected. */
  int check() const ;

protected:
  void refineTriMesh();

  int lookupEdge(const Mesh& mesh, int v1, int v2) const ;
  
private:
  MeshType meshType_;
  Mesh coarse_;
  Mesh fine_;

  Array<int> oldToNewVertMap_;
  Array<int> newVertIsOnEdge_;
  Array<int> newVertToOldLIDMap_;
  Array<int> oldEdgeToNewVertMap_;

  ArrayOfTuples<int> oldToNewElemMap_;
  Array<int> newToOldElemMap_;

  Array<Array<int> > oldEdgeChildren_;
  Array<Array<int> > oldEdgeParallels_;
  Array<int> newEdgeParents_;
  Array<int> newEdgeParallels_;

  ArrayOfTuples<int> interiorEdges_;
};

}


#endif
