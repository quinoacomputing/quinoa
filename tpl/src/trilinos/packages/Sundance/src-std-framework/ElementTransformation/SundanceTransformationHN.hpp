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
 * SundanceTransformationHN.hpp
 *
 *  Created on: Mar 16, 2010
 *      Author: benk
 */

#ifndef SUNDANCETRANSFORMATIONHN_HPP_
#define SUNDANCETRANSFORMATIONHN_HPP_

#include "SundanceTransformationBase.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceNodalDOFMapHN.hpp"
#include "SundanceHNDoFMapBase.hpp"

namespace Sundance {

class TransformationHN: public Sundance::TransformationBase {
public:

	TransformationHN(const HNDoFMapBase* dofMap ,
			const int nrCol , const int nrRaw);

	virtual ~TransformationHN();

	/** */
  virtual void preApply( const int funcID ,
	         int cellDim ,
			 const CellJacobianBatch& JTrans,
			 const CellJacobianBatch& JVol,
			 const Array<int>& facetIndex,
			 const RCP<Array<int> >& cellLIDs,
			 RCP<Array<double> >& A ) const ;

	/** */
  virtual void postApply( const int funcID ,
		      int cellDim ,
			  const CellJacobianBatch& JTrans,
			  const CellJacobianBatch& JVol,
			  const Array<int>& facetIndex,
			  const RCP<Array<int> >& cellLIDs,
			  RCP<Array<double> >& A
			  ) const ;
  
  /** */
  // this will potentially used in scatter process
  virtual void preapplyTranspose( const int cellDim,
				  const int funcID,
				  const Array<int>& cellLIDs,
				  const Array<int>& facetIndex,
				  Array<double>& A ) const;

protected:

	void multiplyFromLeftWithTransp(Array<double>& M , double* A_end , const double* A_copy) const;

	void multiplyFromRight(double* A_end , Array<double>& M , const double* A_copy) const;

	void multiplyFromLeft(Array<double>& M , double* A_end , const double* A_copy ,
			              const int nrRow, const int nrCol) const;

private:

	/** The DoF Map */
	const HNDoFMapBase* dofMap_;

	/** */
	const int nrRow_;

	/** */
	const int nrCol_;

};

}

#endif /* SUNDANCETRANSFORMATIONHN_HPP_ */
