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

#ifndef SUNDANCE_INHOMOGENOUSDOFMAPHN_H
#define SUNDANCE_INHOMOGENOUSDOFMAPHN_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceHNDoFMapBase.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceBasisFamilyBase.hpp"
#include "SundanceMatrixStore.hpp"

namespace Sundance
{
using namespace Teuchos;

class BasisDOFTopologyBase;

/** 
 * A MixedDOFMapHN is a DOF map for the case where
 * every function is defined
 * on every cell in the mesh, but where functions may have different bases. <br>
 * MixedDOFMapHN is defined on meshes with hanging nodes,
 * it builds the restriction matrix where it is needed.
 */
class InhomogeneousDOFMapHN : public HNDoFMapBase , public DOFMapBase
{
public:
  /** */
	InhomogeneousDOFMapHN(const Mesh& mesh,
    const Array<RCP<BasisDOFTopologyBase> >& basis,
    const Array<Map<Set<int>, CellFilter> >& funcSetToDomainMap,
    int setupVerb);
                        
  /** */
  virtual ~InhomogeneousDOFMapHN(){;}

  /** */
  RCP<const MapStructure> 
  getDOFsForCellBatch(int cellDim,
    const Array<int>& cellLID,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verbosity) const ;

	/** Returns the
	 * @param cellLID [in] the maxCell LID input
	 * @param funcID [in]
	 * @param trafoMatrixSize [in/out]
	 * @param doTransform [out]
	 * @param transfMatrix [out] row major matrix*/
  void getTrafoMatrixForCell(
	    int cellLID,
	    int funcID,
	    int& trafoMatrixSize,
	    bool& doTransform,
	    Array<double>& transfMatrix ) const;

	/** Function to apply transformation for facets
	 * @param cellDim , the facet dimension
	 * @param cellLID , facet LID
	 * @param facetIndex , facet index in the maxCofacet
	 * @param funcID  [in] the function ID
	 * @param trafoMatrixSize [in/out]
	 * @param doTransform [out]
	 * @param transfMatrix [out] (we assume that the array is already pre-sized )*/
     void getTrafoMatrixForFacet(
		  int cellDim,
		  int cellLID,
		  int facetIndex,
		  int funcID,
		  int& trafoMatrixSize,
		  bool& doTransform,
		  Array<double>& transfMatrix ) const;

   /** Return the set of a function IDs for a batch of cells for
   * which DOFs can be obtained */
   RCP<const Set<int> >
   allowedFuncsOnCellBatch(int cellDim,
       const Array<int>& cellLID) const ;

  /** See subclass for docu
   * We can be sure that this is used only for nodal data (for VTK plotting)<br>
   * This is implemented only for Nodal data call !!! */
  void getDOFsForHNCell(
	  int cellDim,
	  int cellLID,
      int funcID,
      Array<int>& dofs ,
      Array<double>& coefs ) const;

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

  /** Return the function ID for one chunk and index, same as the previous but returns one value*/
  int getfuncID(int basisChunk, int funcNr ) const{
	  return (this->funcID(basisChunk))[funcNr];
  }

  /** */
  const Array<CellFilter>& funcDomains() const {return funcDomains_;}

  /** */
  virtual void print(std::ostream& os) const ;

private:

  /** */
  inline int getInitialDOFForCell(int cellDim, int cellLID, int funcID) const
    {
      return dofs_[cellDim][cellLID][funcID][0];
    }

  inline int* getInitialDOFPtrForCell(int cellDim, int cellLID, int funcID)
    {
      return &(dofs_[cellDim][cellLID][funcID][0]);
    }

  inline const int* getInitialDOFPtrForCell(int cellDim, int cellLID, 
    int funcID) const
    {
      return &(dofs_[cellDim][cellLID][funcID][0]);
    }

  /** */
  void allocate(const Mesh& mesh);
      
  /** */
  void buildMaximalDofTable();

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

  /** spatial dimension */
  int dim_;

  /** Tables of DOFs, indexed by dimension and chunk number.
   *
   * dof(cellDim, cellLID, funcID , node)
   * = dofs_[cellDim][cellLID][funcID][node]
   */
  Array< Array<Array<Array<int> > > > dofs_;

  /** if one function is defined on this element then is > 0 ,
   * else the function is not defined on this element
   * = funcDefined_[cellDim][funcID][cellLID]*/
  Array< Array<Array<short signed int> > > funcDefined_;

  /** DOFs for maximal cells, indexed by the cell LID
   * this contains all the DoF of one maxDimCell
   * = maximalDofs_[LID][funcID][DoFIndex];
   */
  mutable Array< Array<Array<int> > > maximalDofs_;

  /** whether maximal DOFs have been tabulated */
  mutable bool haveMaximalDofs_;

  /** 
   * localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber]
   */
  Array<Array<Array<Array<Array<int> > > > > localNodePtrs_;

  /** the domains where each cell filter (as it is in the input) is defined on the maxCells */
  Array<CellFilter> maxSubdomains_;

  /** the domains where each function(input) is defined on the maxCells */
  Array<CellFilter> funcDomains_;

  /** set which contains all the functions */
  Set<int> allFuncs_;

  /** here we store all combinations of functions which can appear for maxcells in the input parameter
   * used only in the setup face*/
  Array<Set<int> > elemFuncSetsDomain_;

  /** store all the combination of functions */
  Array<Set<int> > elemFuncSets_;

  /** how many cells has one combination of "elemFuncSets_" */
  Array<int> funcSetCellCount_;

  /** for each maxcell store the index in the "elemFuncSets_" , to store
   * the combination of functions */
  Array<int> maxCellFuncSetsIndex_;

/* ========= Hanging node treatment ============= */

   /** Nr of points, needed only by plotting */
   int nPoints_;

   /** total number of functions */
   int nrAllFuncs_;

   /** Is true if the cell has hanging node */
   Array<bool> hasCellHanging_;

   /** Flag per element showing if the element is a hanging element*/
   Array< Array<bool> > isElementHanging_;

   /** store the facetLIDs of the points [nPoints*b + facetID]*/
   Sundance::Map< int , Array<int>  > HN_To_globalFacetsLID_;

   /** store the facetDim of the points [nPoints*b + facetID]*/
   Sundance::Map< int , Array<int>  > HN_To_globalFacetsDim_;

   /** store the coefs [nPoints*b + facetID]*/
   Sundance::Map< int , Array<double>  > HN_To_coeffs_;

   /** Maps one maxCell LID to one transformation matrix, only maxDim elements and basisChunck ID
    * with hanging nodes should have transformation matrix
    * [maxCellLID].get -> [basisChunckID][] is the index of trafo matrix */
   Sundance::Map< int , Array < int > > maxCellLIDwithHN_to_TrafoMatrix_;

   /** The object to store all the transformation matrixes */
   MatrixStore                          matrixStore_;

   /** The basis for the transformation */
   mutable Array<RCP<BasisFamilyBase> > basis_;

/* ========= END hanging node treatment ======== */

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
