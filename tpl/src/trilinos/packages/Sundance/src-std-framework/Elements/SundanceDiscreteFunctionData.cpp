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

#include "SundanceDiscreteFunctionData.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceSubtypeEvaluator.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#endif


using namespace Sundance;
using namespace Teuchos;




DiscreteFunctionData::DiscreteFunctionData(const DiscreteSpace& space)
  : DiscreteFuncDataStub(), 
    space_(space),
    vector_(space_.createVector()),
    ghostView_(),
    ghostsAreValid_(false)
{}

DiscreteFunctionData::DiscreteFunctionData(const DiscreteSpace& space, 
  const double& constantValue)
  : DiscreteFuncDataStub(), 
    space_(space),
    vector_(space_.createVector()),
    ghostView_(),
    ghostsAreValid_(false)
{
  vector_.setToConstant(constantValue);
}

DiscreteFunctionData::DiscreteFunctionData(const DiscreteSpace& space, 
  const Vector<double>& vector)
  : DiscreteFuncDataStub(), 
    space_(space),
    vector_(vector),
    ghostView_(),
    ghostsAreValid_(false)
{}

const DiscreteFunctionData* DiscreteFunctionData::getData(const DiscreteFuncElement* dfe)
{
  TEUCHOS_TEST_FOR_EXCEPTION(dfe==0, std::runtime_error, "null argument to DiscreteFunctionData::getData()");
  RCP<const DiscreteFunctionData> rtn 
    = rcp_dynamic_cast<const DiscreteFunctionData>(dfe->commonData());
  TEUCHOS_TEST_FOR_EXCEPTION(rtn.get()==0, std::runtime_error, 
    "cast to DiscreteFunctionData* failed for "
    "discrete function element " << dfe->toXML());
  return rtn.get();
}

DiscreteFunctionData* DiscreteFunctionData::getData(DiscreteFuncElement* dfe)
{
  TEUCHOS_TEST_FOR_EXCEPTION(dfe==0, std::runtime_error, "null argument to DiscreteFunctionData::getData()");
  DiscreteFunctionData* rtn 
    = dynamic_cast<DiscreteFunctionData*>(dfe->commonData());
  TEUCHOS_TEST_FOR_EXCEPTION(rtn==0, std::runtime_error, 
    "cast to DiscreteFunctionData* failed for "
    "discrete function element " << dfe->toXML());
  return rtn;
}

void DiscreteFunctionData::setVector(const Vector<double>& vec) 
{
  ghostsAreValid_ = false;
  vector_ = vec;
}

void DiscreteFunctionData::updateGhosts() const
{
  if (!ghostsAreValid_)
  {
    space_.importGhosts(vector_, ghostView_);
    ghostsAreValid_ = true;
  }
}


RCP<const MapStructure> DiscreteFunctionData
::getLocalValues(int cellDim, 
  const Array<int>& cellLID,
  Array<Array<double> >& localValues) const 
{
  Tabs tab;

  if (Evaluator::classVerbosity() > 3)
  {
    Out::os() << tab << "getting DF local values" << std::endl;
  }
  updateGhosts();

  const RCP<DOFMapBase>& map = space_.map();
  Array<Array<int> > dofs;
  Array<int> nNodes;

  RCP<const Set<int> > requestedFuncs = map->allowedFuncsOnCellBatch(cellDim,
    cellLID);

  RCP<const MapStructure> s = map->getDOFsForCellBatch(cellDim, cellLID,
    *requestedFuncs,
    dofs, nNodes, Evaluator::classVerbosity());
  localValues.resize(s->numBasisChunks());
  // test if we need any kind of transformation

  if  (space_.getTransformation()->validTransformation())
  {
	 SUNDANCE_OUT( Evaluator::classVerbosity() > 3 , "DiscreteFunctionData::getLocalValues() VALID TRAFO FOUND ... ")
	 for (int b=0; b<nNodes.size(); b++)
	 {
	    int nFuncs = s->numFuncs(b);
	    Array<int> functionIDs = s->funcs(b);
	    localValues[b].resize(nFuncs*nNodes[b]*cellLID.size());

	    // first get the dofs values, which later will be transformed
        ghostView_->getElements(&(dofs[b][0]), dofs[b].size(), localValues[b]);
        // make transformation and fill the correct "localValues[b]" elements !!!! (if it is needed)
	    // do the transformation for each function , ("nFuncs")
        // nNodes[b] is the total number for one function inside the chunk
	    space_.getTransformation()->getDoFsWithTransformation(
	    		dofs[b] , functionIDs , b , nNodes[b] , nFuncs , cellDim, cellLID ,ghostView_ , localValues[b] );
	 }
  }
  else  // if we do not need transformation then do the normal thing
  {
	  SUNDANCE_OUT( Evaluator::classVerbosity() > 3 , "DiscreteFunctionData::getLocalValues() NO VALID TRAFO ... ")
	  for (int b=0; b<nNodes.size(); b++)
	  {
	    int nFuncs = s->numFuncs(b);
	    localValues[b].resize(nFuncs*nNodes[b]*cellLID.size());
        ghostView_->getElements(&(dofs[b][0]), dofs[b].size(), localValues[b]);
	  }
  }

  return s;
}



