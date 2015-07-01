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


/*
 * SundanceAssemblyTransformationBuilder.cpp
 *
 *  Created on: Mar 16, 2010
 *      Author: benk
 */

#include "SundanceAssemblyTransformationBuilder.hpp"
#include "SundanceHNDoFMapBase.hpp"
#include "SundanceMixedDOFMap.hpp"

using namespace Sundance;

/* for Assembly transformations*/
AssemblyTransformationBuilder::AssemblyTransformationBuilder( const RCP<IntegralGroup>& group ,
							      const Array<RCP<DOFMapBase> >& rowMaps ,
							      const Array<RCP<DOFMapBase> >& colMaps ,
							      const Mesh& mesh):
  verb_(0),
  nrCol_(group->nUnkNodes()) ,
  nrRow_(group->nTestNodes()) ,
  preTransformation_(0),
  postTransformation_(0),
  testFuncID_(0),
  unkFuncID_(0),
  _myRowDOFMap(0),
  _myColDOFMap(0),
  hasTransformation_(false),
  hasPreTransformation_(false),
  hasPostTransformation_(false)
{
  SUNDANCE_MSG2(verb(), "AssemblyTransformationBuilder::AssemblyTransformationBuilder initialized fields:" << mesh.allowsHangingHodes());
  // make different transformations
  if (mesh.allowsHangingHodes()){

    // we have HANGING NODES, create a corresponding transformation
    hasTransformation_ = true;
    if (group->nUnkNodes() > 0) // 2-form
      {
	SUNDANCE_MSG2(verb(), "Matrix trafo colMaps.size:" << colMaps.size());
	SUNDANCE_MSG2(verb(), "Matrix trafo rowMaps.size:" << rowMaps.size());
	SUNDANCE_MSG2(verb(), "Matrix trafo group->unkBlock():" << group->unkBlock());
	SUNDANCE_MSG2(verb(), "Matrix trafo group->testBlock():"<< group->testBlock());
	// this means we have to multiply Matrix
	hasPreTransformation_ = true;
	hasPostTransformation_ = true;
	_myColDOFMap = colMaps[group->unkBlock()[0]].get(); // group->unkBlock(), has always one element ;-)
	_myRowDOFMap = rowMaps[group->testBlock()[0]].get(); // group->testBlock(), has always one element ;-)
	SUNDANCE_MSG2(verb(), "Trafo the Maps to HN Maps" );
	// create the NH transformation, for pre and post apply
	const HNDoFMapBase* colMap
	  = dynamic_cast<const HNDoFMapBase*>(_myColDOFMap);
	const HNDoFMapBase* rowMap
	  = dynamic_cast<const HNDoFMapBase*>(_myRowDOFMap);

	testFuncID_ = group->testID()[0];
	unkFuncID_ = group->unkID()[0];

	if (colMap != 0)
	  {
	    postTransformation_ = rcp((TransformationBase*)(new TransformationHN( colMap , nrCol_ , nrRow_ )));
	  }
	else
	  {
	    SUNDANCE_ERROR("AssemblyTransformationBuilder::AssemblyTransformationBuilder, wrong configuration , no Col DoFMap");
	  }
	
	if (rowMap != 0)
	  {
	    preTransformation_ = rcp(new TransformationHN( rowMap , nrCol_ , nrRow_ ));
	  }
	else
	  {
	    SUNDANCE_ERROR("AssemblyTransformationBuilder::AssemblyTransformationBuilder, wrong configuration , no Row DoFMap");
	  }
      }
    else
      {
	// this means we have to multiply Vector
	if (group->nTestNodes() > 0){
	  testFuncID_ = group->testID()[0];
	  SUNDANCE_MSG2(verb(), "Vector trafo colMaps.size:" << colMaps.size());
	  SUNDANCE_MSG2(verb(), "Vector trafo rowMaps.size:" << rowMaps.size());
	  SUNDANCE_MSG2(verb(), "Vector trafo group->unkBlock():" << group->unkBlock());
	  SUNDANCE_MSG2(verb(), "Vector trafo group->testBlock():"<< group->testBlock());
	  _myRowDOFMap = rowMaps[group->testBlock()[0]].get(); // group->testBlock(), has always one element ;-)
	  hasPreTransformation_ = true;
	  hasPostTransformation_ = false;
	  // create the NH transformation, only for pre apply
	  const HNDoFMapBase* rowMap
	    = dynamic_cast<const HNDoFMapBase*>(_myRowDOFMap);

	  if (rowMap != 0)
	    {
	      preTransformation_ = rcp(new TransformationHN( rowMap , nrCol_ , nrRow_ ));
	    }
	  else
	    {
	      SUNDANCE_ERROR("AssemblyTransformationBuilder::AssemblyTransformationBuilder, wrong configuration , no Col DoFMap ");
	    }
	}
	else
	  {
	    hasTransformation_ = false;
	  }
      }
  }
  else
    {
      if (group->nUnkNodes() > 0)  // two-form
	{
	  _myColDOFMap = colMaps[group->unkBlock()[0]].get(); // group->unkBlock(), has always one element ;-)
	  _myRowDOFMap = rowMaps[group->testBlock()[0]].get(); // group->testBlock(), has always one element ;-)
	  const MixedDOFMap* colMap
	    = dynamic_cast<const MixedDOFMap*>(_myColDOFMap);
	  const MixedDOFMap* rowMap
	    = dynamic_cast<const MixedDOFMap*>(_myRowDOFMap);

	  if (colMap != 0) 
	    {
	      postTransformation_ = rcp( new InequivalentElementTransformation( mesh , colMap ) );
	      hasPostTransformation_ = true;
	      hasTransformation_ = true;
	    }
	  if (rowMap != 0)
	    {
	      preTransformation_ = rcp( new InequivalentElementTransformation( mesh , rowMap ) );
	      hasPreTransformation_ = true;
	      hasTransformation_ = true;
	    }
	}
      else if (group->nTestNodes() > 0) // one-form
	{
	  _myRowDOFMap = rowMaps[group->testBlock()[0]].get(); // group->testBlock(), has always one element ;-)
	  const MixedDOFMap* rowMap
	    = dynamic_cast<const MixedDOFMap*>(_myRowDOFMap);

	  if (rowMap != 0) 
	    {
	      preTransformation_ = rcp( new InequivalentElementTransformation( mesh , rowMap ) );
	      hasPreTransformation_ = true;
	      hasTransformation_ = true;
	    }
	}

    }
}

AssemblyTransformationBuilder::~AssemblyTransformationBuilder() {

}

void AssemblyTransformationBuilder::applyTransformsToAssembly( int groupIndex ,
							       int entryPerCell ,
							       CellType cellType ,
							       int cellDim ,
							       CellType maxCellType ,
							       const CellJacobianBatch& JTrans,
							       const CellJacobianBatch& JVol,
							       const Array<int>& facetNum,
							       const RCP<Array<int> >& cellLIDs,
							       RCP<Array<double> >& A)
{
  
  // if the integration is not on a MaxCell then no HangingNode Trafo!
  // caz in some cases this would incease the size of the elem matrix
  //if (cellType != maxCellType){
    //hasTransformation_ = false;
    // in 3D for boundary faces we need a 2D transformation
    // eventhrough boundary faces might be hanging faces, we might need to integrate on them !!!
    //return;
  //}

  if (hasTransformation_)
    {
      if (hasPreTransformation_) 
	{
	  preTransformation_->preApply( testFuncID_ , cellDim ,
					JTrans, JVol, facetNum, cellLIDs, A);
	}
      if (hasPostTransformation_)
	{
	  postTransformation_->postApply( unkFuncID_ , cellDim ,
					  JTrans, JVol, facetNum, cellLIDs, A);
	}
    }
  
  SUNDANCE_MSG2( verb() , " AssemblyTransformationBuilder::applyTransformsToAssembly ");
}
