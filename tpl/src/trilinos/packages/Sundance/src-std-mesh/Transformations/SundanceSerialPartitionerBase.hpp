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

#ifndef SUNDANCE_SERIALPARTITIONERBASE_H
#define SUNDANCE_SERIALPARTITIONERBASE_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshType.hpp"
#include "PlayaHandle.hpp"

namespace Sundance
{


/**
 * Base class for mesh partitioners that run in serial
 */
class SerialPartitionerBase
{
public:
  /** */
  SerialPartitionerBase(bool ignoreGhosts=false)
    : ignoreGhosts_(ignoreGhosts){}
  
  /** */
  virtual ~SerialPartitionerBase(){;}

  /** */
  void getNeighbors(const Mesh& mesh, 
    Array<Array<int> >& neighbors, int& nEdges) const ;

  /** */
  Set<int> arrayToSet(const Array<int>& a) const ;

  /** */
  virtual void getAssignments(const Mesh& mesh, int np, 
    Array<int>& assignments) const = 0 ;

  /** */
  Array<Mesh> makeMeshParts(const Mesh& mesh, int np,
    Array<Sundance::Map<int, int> >& oldElemLIDToNewLIDMap,
    Array<Sundance::Map<int, int> >& oldVertLIDToNewLIDMap
    ) const ;

  /** */
  void getOffProcData(int p, 
    const Array<int>& elemAssignments,
    const Array<int>& nodeAssignments,
    Set<int>& offProcNodes,
    Set<int>& offProcElems) const ;

  /** 
   * 
   */
  void getNodeAssignments(int nProc, 
    const Array<int>& elemAssignments,
    Array<int>& nodeAssignments,
    Array<int>& nodeOwnerElems,
    Array<int>& nodesPerProc) const ;

  /** */
  void getElemsPerProc(int nProc, 
    const Array<int>& elemAssignments,
    Array<int>& elemsPerProc) const ;

  /** Remap global element or node 
   * numberings so that each processor owns sequentially-numbered
   * global indexes. */
  void remapEntities(const Array<int>& assignments, int nProc,
    Array<int>& entityMap) const ;


private:

  bool ignoreGhosts_;
  int max(const Set<int>& s) const ;
  mutable Array<Set<int> > elemVerts_;
  mutable Array<Set<int> > elemEdgewiseNbors_;
  mutable Array<Set<int> > vertElems_;
};
}

#endif
