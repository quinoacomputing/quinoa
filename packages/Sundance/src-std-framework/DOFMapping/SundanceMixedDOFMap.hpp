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

#ifndef SUNDANCE_MIXEDDOFMAP_H
#define SUNDANCE_MIXEDDOFMAP_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceSpatiallyHomogeneousDOFMapBase.hpp"

namespace Sundance
{
using namespace Teuchos;

class BasisDOFTopologyBase;

/** 
 * A MixedDOFMap is a DOF map for the case where 
 * every function is defined
 * on every cell in the mesh, but where functions may have different bases. 
 */
class MixedDOFMap : public SpatiallyHomogeneousDOFMapBase
{
public:
  /** */
  MixedDOFMap(const Mesh& mesh, 
    const Array<RCP<BasisDOFTopologyBase> >& basis,
    const CellFilter& maxCells, 
    int setupVerb);
                        
  /** */
  virtual ~MixedDOFMap(){;}

  /** */
  RCP<const MapStructure> 
  getDOFsForCellBatch(int cellDim,
    const Array<int>& cellLID,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verbosity) const ;

  /** */
  RCP<const MapStructure> mapStruct() const 
    {return structure_;}

  /** */
  int chunkForFuncID(int funcID) const
    {return structure_->chunkForFuncID(funcID);}

  /** */
  int indexForFuncID(int funcID) const 
    {return structure_->indexForFuncID(funcID);}
      
  /** */
  int nFuncs(int basisChunk) const
    {return nFuncs_[basisChunk];}

  /** */
  int nBasisChunks() const 
    {return nFuncs_.size();}

  /** */
  const RCP<BasisDOFTopologyBase>& basis(int basisChunk) const
    {return structure_->basis(basisChunk);}

  /** */
  const Array<int>& funcID(int basisChunk) const 
    {return structure_->funcs(basisChunk);}


private:

  /** */
  void checkTable() const ;

  /** */
  inline int getInitialDOFForCell(int cellDim, int cellLID, int basisChunk) const
    {
      return dofs_[cellDim][basisChunk][cellLID*nDofsPerCell_[basisChunk][cellDim]];
    }

  inline int* getInitialDOFPtrForCell(int cellDim, int cellLID, int basisChunk)
    {
      return &(dofs_[cellDim][basisChunk][cellLID*nDofsPerCell_[basisChunk][cellDim]]);
    }

  inline const int* getInitialDOFPtrForCell(int cellDim, int cellLID, 
    int basisChunk) const 
    {
      return &(dofs_[cellDim][basisChunk][cellLID*nDofsPerCell_[basisChunk][cellDim]]);
    }

  /** */
  void allocate(const Mesh& mesh);
      
  /** */
  void buildMaximalDofTable() const ;

  bool hasBeenAssigned(int cellDim, int cellLID) const 
    {return hasBeenAssigned_[cellDim][cellLID];}

  void markAsAssigned(int cellDim, int cellLID)
    {hasBeenAssigned_[cellDim][cellLID] = true;}

  /** */
  void initMap();

  /** */
  void setDOFs(int basisChunk, int cellDim, int cellLID, 
    int& nextDOF, bool isRemote=false);

  /** */
  void shareDOFs(int cellDim,
    const Array<Array<int> >& outgoingCellRequests);

  /** */
  void computeOffsets(int dim, int localCount);

  /** */
  static int uninitializedVal() {return -1;}

  /** */
  CellFilter maxCells_;

  /** spatial dimension */
  int dim_;

  /** Tables of DOFs, indexed by dimension and chunk number.
   *
   * dof(cellDim, cellLID, chunk, func, node) 
   * = dofs_[cellDim][chunk][(cellLID*nFunc + func)*nNode + node]
   */
  Array<Array<Array<int> > > dofs_;

  /** DOFs for maximal cells, indexed by basis chunk number 
   *
   * dof(cellLID, chunk, func, node) 
   * = maximalDofs_[chunk][(cellLID*nFunc + func)*nNode + node];
   */
  mutable Array<Array<int> > maximalDofs_;

  /** whether maximal DOFs have been tabulated */
  mutable bool haveMaximalDofs_;

  /** 
   * localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber]
   */
  Array<Array<Array<Array<Array<int> > > > > localNodePtrs_;

  /** The number of nodes per cell, for each basis function type, 
   * <i>not</i> including the nodes of the facets of the cell. Indexed as 
   * nNodesPerCell_[basis][dimension] */
  Array<Array<int> > nNodesPerCell_;

  /** The number of DOFs per cell, for each basis function type, 
   * <i>not</i> including the DOFs of the facets of the cell. Indexed as 
   * nDofsPerCell_[basis][dimension] */
  Array<Array<int> > nDofsPerCell_;

  /** The number of nodes per cell, for each basis function type, including
   * the nodes of the facets of the cell. Indexed as 
   * nNodesPerCell_[basis][dimension] */
  Array<Array<int> > totalNNodesPerCell_;

  /** The number of DOFs per cell, for each basis function type, including
   * the DOFs of the facets of the cell. Indexed as 
   * nDofsPerCell_[basis][dimension] */
  Array<Array<int> > totalNDofsPerCell_;

  /** Indicates whether the cells of each dimension have any DOFs in this
   * map, for any chunk. */
  Array<int> cellHasAnyDOFs_;

  /** number of facets of dimension facetDim for cells of dimension cellDim.
   * Indexed as numFacets_[cellDim][facetDim]
   */
  Array<Array<int> > numFacets_;

  /** Orientation of each edge or face as seen by the maximal cell
   * from which its DOFs were originally assigned. */
  Array<Array<int> > originalFacetOrientation_;

  /** */
  Array<Array<int> > hasBeenAssigned_;

  /** */
  RCP<const MapStructure> structure_;

  /** */
  Array<int> nFuncs_;
};
}

                  
#endif
