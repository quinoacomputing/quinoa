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
 * SundanceInequivalentTransformation.hpp
 *
 *  Created on: 5 April, 2010
 *      Author: R.C. Kirby
 */

#ifndef SUNDANCEINEQUIVALENTELEMENTTRANSFORMATION_HPP_
#define SUNDANCEINEQUIVALENTELEMENTTRANSFORMATION_HPP_

#include "SundanceIntegralGroup.hpp"
#include "SundanceTransformationBase.hpp"
#include "SundanceMixedDOFMap.hpp"

namespace Sundance {

  class InequivalentElementTransformation : public TransformationBase {
  public:
    
    InequivalentElementTransformation( const Mesh &mesh, 
				       const MixedDOFMap *map );
    
    virtual ~InequivalentElementTransformation() {;}
    
    /** The transformation method */
    // this will potentially used in assembly  process
    virtual void preApply( const int funcID,
    		   int cellDim ,
			   const CellJacobianBatch& JTrans,
			   const CellJacobianBatch& JVol,
			   const Array<int>& facetIndex,
			   const RCP<Array<int> >& cellLIDs,
			   RCP<Array<double> >& A
			   ) const;
    
    /** */
    // this will potentially used in assembly  process
    virtual void postApply( const int funcID,
    		    int cellDim ,
			    const CellJacobianBatch& JTrans,
			    const CellJacobianBatch& JVol,
			    const Array<int>& facetIndex,
			    const RCP<Array<int> >& cellLIDs,
			    RCP<Array<double> >& A
			    ) const;

    /** */
    // this will potentially used in scatter process
    virtual void preapplyTranspose( const int cellDim,
				    const int funcID,
				    const Array<int>& cellLIDs,
				    const Array<int>& facetIndex,
				    Array<double>& A ) const;



  protected:

    /** */
    int verb() const { return verb_; }

    /** */
    void setverb(int c) { verb_ = c; }

  private :
    /** reference to the mesh, needed for preapply transpose interface */
    const Mesh &mesh_;

    /** reference to the basis; transformation calls are forwarded to it */
    const MixedDOFMap *map_;

    /** pointers to bases for each chunk */
    Array<const BasisFamilyBase *> chunkBases_;

    /** verbosity attribute */
    int verb_;
  };

}

#endif /* SUNDANCEINEQUIVALENTELEMENTTRANSFORMATION_HPP_ */
