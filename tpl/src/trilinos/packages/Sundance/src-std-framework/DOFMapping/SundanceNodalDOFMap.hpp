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

#ifndef SUNDANCE_NODALDOFMAP_H
#define SUNDANCE_NODALDOFMAP_H

#include "SundanceDefs.hpp"
#include "SundanceSpatiallyHomogeneousDOFMapBase.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceObjectWithVerbosity.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * 
 */
class NodalDOFMap : public SpatiallyHomogeneousDOFMapBase
{
public:
  /** */
  NodalDOFMap(const Mesh& mesh, int nFuncs,
    const CellFilter& maxCellFilter,
    int setupVerb);
      
  /** */
  virtual ~NodalDOFMap(){;}

  /** */
  RCP<const MapStructure> 
  getDOFsForCellBatch(int cellDim,
    const Array<int>& cellLID,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verb) const ;

  /** */
  RCP<const MapStructure> mapStruct() const 
    {return structure_;}

  /** */
  int nFuncs() const {return nFuncs_;}

protected:

  void init();

  void computeOffsets(int localCount)  ;

  void shareRemoteDOFs(const Array<Array<int> >& remoteNodes);

  CellFilter maxCellFilter_;
  int dim_;
  int nFuncs_;
  int nElems_;
  int nNodes_;
  int nNodesPerElem_;
  Array<int> elemDofs_;
  Array<int> nodeDofs_;
  RCP<const MapStructure> structure_;
};

}


#endif
