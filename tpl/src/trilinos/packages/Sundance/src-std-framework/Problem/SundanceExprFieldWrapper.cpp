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

#include "SundanceExprFieldWrapper.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceEdgeLocalizedBasis.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceHNDoFMapBase.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace Playa;


ExprFieldWrapper::ExprFieldWrapper(const Expr& expr)
  : expr_(expr),
    df_(),
    discreteSpace_(),
    //map_(),
    indices_(),
    Expr_size_(1),
    isPointData_(true)
{
	int index = 0;
	Expr_size_ = expr.size();
	// Now it is independent of the size of the size of the Expression
	for(index = 0 ; index < Expr_size_ ; index++)
	{
	  const DiscreteFunction* df
      = dynamic_cast<const DiscreteFunction*>(expr[index].ptr().get());
	  const DiscreteFuncElement* dfe
      = dynamic_cast<const DiscreteFuncElement*>(expr[index].ptr().get());
    if (df != 0)
    {
      discreteSpace_ = df->discreteSpace();
      //map_ = df->map();
      indices_.append(tuple(0));
      BasisFamily basis = discreteSpace_.basis()[0];
      const Lagrange* lagr = dynamic_cast<const Lagrange*>(basis.ptr().get());
      if (lagr != 0 && lagr->order()==0) isPointData_ = false;
      const EdgeLocalizedBasis* elb = dynamic_cast<const EdgeLocalizedBasis*>(basis.ptr().get());
      if (elb!=0) isPointData_ = false;
      df_ = df->data();
    }
    else if (dfe != 0)
    {
      const DiscreteFunctionData* f = DiscreteFunctionData::getData(dfe);

      TEUCHOS_TEST_FOR_EXCEPTION(f == 0, std::runtime_error,
        "ExprFieldWrapper ctor argument "
        << expr << " is not a discrete function");
      discreteSpace_ = f->discreteSpace();
      //map_ = f->map();
      indices_.append(tuple(dfe->myIndex()));
      BasisFamily basis = discreteSpace_.basis()[indices_[index][0]];
      const Lagrange* lagr = dynamic_cast<const Lagrange*>(basis.ptr().get());
      if (lagr != 0 && lagr->order()==0) isPointData_ = false;      
      const EdgeLocalizedBasis* elb = dynamic_cast<const EdgeLocalizedBasis*>(basis.ptr().get());
      if (elb!=0) isPointData_ = false;

      df_ = f;
          
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(df == 0 && dfe == 0, std::runtime_error,
        "ExprFieldWrapper ctor argument is not a discrete "
        "function");
    }
  }
}


double ExprFieldWrapper::getData(int cellDim, int cellID, int elem) const
{
  Array<int> dofs;

  discreteSpace_.map()->getDOFsForCell(cellDim, cellID, indices_[elem][0] , dofs); //indecies[elem][0] should be OK!

  //cout << "Arguments ExprFieldWrapper::getData " << cellDim << "," << cellID << "," << elem << " DoFs:" << dofs <<std::endl;
  //cout << "indices_[elem][0]:" << indices_[elem][0] << std::endl;

  // This exception is not needed since the first value must be the nodal value
  //TEUCHOS_TEST_FOR_EXCEPTION(dofs.size() > 1, std::runtime_error,
  //  "too many DOFs found in ExprFieldWrapper::getData()");

  // in case of hanging node the "dofs[0]" will be negative, in this case treate this case
  // in case of general basis function this should not be changed, when there are nodal values
  // Todo: if we do not have nodal values then we should do some extra things ...
  if ( dofs[0] < 0)
  {
  	const HNDoFMapBase* HNMap
  		    = dynamic_cast<const HNDoFMapBase*>((discreteSpace_.map()).get());
    if (HNMap != 0 ){
        Array<double> coefs;
        double sum = 0.0;
    	HNMap->getDOFsForHNCell( cellDim, cellID, indices_[elem][0] ,  dofs , coefs );
    	for (int jj = 0 ; jj < dofs.size() ; jj++)
    	{
    		sum += coefs[jj] * df_->ghostView()->getElement(dofs[jj]); //sum up the contributions
    	}
    	// return the contribution of the global DoFs
    	return sum;
    	//return 1.0;
    }
    else
    {
	  return 1.0;
    }
  }
  else
  {
	  return df_->ghostView()->getElement(dofs[0]);
  }
}


bool ExprFieldWrapper::isDefined(int cellDim, int cellID, int elem) const
{
  // this works only for the first
  RCP<const Set<int> > allowedFuncs 
    = discreteSpace_.map()->allowedFuncsOnCellBatch(cellDim, tuple(cellID));

  //cout << "Arguments ExprFieldWrapper::isDefined" << cellDim << "," << cellID << "," << elem << std::endl;

  return allowedFuncs->contains(indices_[elem][0]);
}
