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


#include "SundanceInequivalentElementTransformation.hpp"

namespace Sundance
{
  InequivalentElementTransformation::InequivalentElementTransformation( const Mesh& mesh ,
									const MixedDOFMap *map ):
    mesh_( mesh ), map_( map ), chunkBases_( map->nBasisChunks() )
  {
	bool isAnyTransfReq = false;
    // extract all the bases from the dof map
    for (int i=0;i<map->nBasisChunks();i++) 
      {
	chunkBases_[i] = dynamic_cast<const BasisFamilyBase *>(map->basis(i).get());

	if (chunkBases_[i] != 0) { isAnyTransfReq = (isAnyTransfReq || chunkBases_[i]->requiresBasisTransformation()); }
      }
    // set the flag which shows if there is any transformation needed in the array of basis
    setDoesAnyTransformation(isAnyTransfReq);
  }

  void InequivalentElementTransformation::preApply( const int funcID,
		                    int cellDim ,
						    const CellJacobianBatch& JTrans,
						    const CellJacobianBatch& JVol,
						    const Array<int>& facetIndex,
						    const RCP<Array<int> >& cellLIDs,
						    RCP<Array<double> >& A
						    ) const
  
  {
	// do integration for MaxDimCells
	if ( mesh_.spatialDim() > cellDim) return;

	if (chunkBases_[map_->chunkForFuncID( funcID )]->requiresBasisTransformation())
	{
	   CellJacobianBatch JVol1;
	   mesh_.getJacobians( mesh_.spatialDim() , *(cellLIDs.get()) , JVol1 );
	   chunkBases_[map_->chunkForFuncID( funcID )]->preApplyTransformation( mesh_.cellType( mesh_.spatialDim() ) , mesh_, *cellLIDs, JVol1 , A );
	}
  }

  void InequivalentElementTransformation::postApply( const int funcID,
		                     int cellDim ,
						     const CellJacobianBatch& JTrans,
						     const CellJacobianBatch& JVol,
						     const Array<int>& facetIndex,
						     const RCP<Array<int> >& cellLIDs,
						     RCP<Array<double> >& A
						     ) const
  {
	// do integration for MaxDimCells
	if ( mesh_.spatialDim() > cellDim) return;

	if (chunkBases_[map_->chunkForFuncID( funcID )]->requiresBasisTransformation())
	{
	   CellJacobianBatch JVol1;
	   mesh_.getJacobians( mesh_.spatialDim() , *(cellLIDs.get()) , JVol1 );
	   chunkBases_[map_->chunkForFuncID( funcID )]->postApplyTransformation( mesh_.cellType( mesh_.spatialDim() ) , mesh_ , *cellLIDs, JVol1 , A );
	}
  }
  
  void InequivalentElementTransformation::preapplyTranspose( const int cellDim,
							     const int funcID,
							     const Array<int>& cellLIDs,
							     const Array<int>& facetIndex,
							     Array<double>& A ) const
  {
	// do the transformation only when it is needed so we do not have to query the Jacobians
	if (chunkBases_[map_->chunkForFuncID( funcID )]->requiresBasisTransformation())
	{
		CellJacobianBatch JVol;
		mesh_.getJacobians( mesh_.spatialDim() , cellLIDs , JVol );
		chunkBases_[map_->chunkForFuncID( funcID )]->preApplyTransformationTranspose( mesh_.cellType( mesh_.spatialDim() ) , mesh_ ,  cellLIDs, JVol , A );
	}
  }
  

}
