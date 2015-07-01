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
 * SundanceDiscreteSpaceTransfBuilder.cpp
 *
 *  Created on: Mar 21, 2010
 *      Author: benk
 */


#include "SundanceDiscreteSpaceTransfBuilder.hpp"


using namespace Sundance;


DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder():
  verb_(0),
  basisSize_(0),
  hasTransformation_(false),
  transformation_()
{
}

DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder( const Mesh& mesh, const BasisArray& basis,
							const RCP<DOFMapBase>& map):
  verb_(0),
  basisSize_(0),
  hasTransformation_(false),
  transformation_()
{

  if (mesh.allowsHangingHodes())
    {
      // in this case we have to create a transformation
      const HNDoFMapBase* myMap
	= dynamic_cast<const HNDoFMapBase*>(map.get());
      if (myMap != 0)
	{
	  // create one hanging node transformation
	  SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder basis.size():" << basis.size());
	  basisSize_ = basis.size();
	  transformation_ = rcp((TransformationBase*)(new TransformationHN( myMap , 1 , basis.size() ))); //todo: basis.size is wrong
	  // the transformation is always needed when we have a mesh with hanging nodes
	  hasTransformation_ = true;
	}
    }
  else
    {
      const MixedDOFMap * myMap = dynamic_cast<const MixedDOFMap *>( map.get() );
      if (myMap != 0)
	{
	  SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::DiscreteSpaceTransfBuilder basis.size():" << basis.size());
	  basisSize_ = basis.size();
	  transformation_ = rcp((TransformationBase*)(new InequivalentElementTransformation( mesh , myMap )));
	  hasTransformation_ = transformation_->doesAnyTransformation();
	}
    }
}

void DiscreteSpaceTransfBuilder::getDoFsWithTransformation( const Array<int>& dofs,
							    const Array<int>& functionIDs ,
							    const int chunkID ,
							    const int nrDoFsPerCell ,
							    const int nrFunctions ,
							    const int cellDim ,
							    const Array<int>& cellLIDs,
							    const RCP<GhostView<double> >& ghostView ,
							    Array<double>& localValues) const 
{
  if (hasTransformation_) 
    {
      SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation()" );
      // get all the DoFs for this functionID
      int DoFPerElement = nrDoFsPerCell;
      Array<double> tmpArray(DoFPerElement*cellLIDs.size());
      Array<int> facetIndex(1); //just needed for the interface
      int cellI , elemDof;
      
      SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation()  nrFunctions:" << nrFunctions
		     << " DoFPerElement:" << DoFPerElement <<  "  localValues.size():" << localValues.size());
      SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation() nrDoFsPerCell:" << nrDoFsPerCell);
      SUNDANCE_MSG2( verb() , "DiscreteSpaceTransfBuilder::getDoFsWithTransformation() localValues:" << localValues);
      
      for (int nf = 0 ; nf < nrFunctions ; nf++)
	{
	  // copy the element values into the temporary array
	  for(cellI = 0 ; cellI < cellLIDs.size() ; cellI++) {
	    for (elemDof = 0 ; elemDof < DoFPerElement ; elemDof++)
	      tmpArray[cellI*DoFPerElement + elemDof] = localValues[(nrFunctions*cellI+nf)*DoFPerElement + elemDof];
	  }
	   SUNDANCE_MSG2( verb() ,"getDoFsWithTransformation() before Transformation:" << tmpArray);
	   SUNDANCE_MSG2( verb() ,"getDoFsWithTransformation() chunk:" << chunkID <<" functionID:" << functionIDs[nf]);
      // make the transformation for all elements once
	   transformation_->preapplyTranspose(
			 cellDim ,
			 functionIDs[nf] ,
			 cellLIDs ,
			 facetIndex,
			 tmpArray );
	   // copy the element values back
	   SUNDANCE_MSG2( verb() , "getDoFsWithTransformation() after Transformation:" << tmpArray );

	  for(cellI = 0 ; cellI < cellLIDs.size() ; cellI++) {
	    for (elemDof = 0 ; elemDof < DoFPerElement ; elemDof++)
	      localValues[(nrFunctions*cellI+nf)*DoFPerElement + elemDof] = tmpArray[cellI*DoFPerElement + elemDof];
	  }
	}
    }
}
