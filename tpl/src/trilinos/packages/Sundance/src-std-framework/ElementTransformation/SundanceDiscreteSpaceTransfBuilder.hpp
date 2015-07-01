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
 * SundanceDiscreteSpaceTransfBuilder.hpp
 *
 *  Created on: Mar 21, 2010
 *      Author: benk
 */

#ifndef SUNDANCEDISCRETESPACETRANSFBUILDER_HPP_
#define SUNDANCEDISCRETESPACETRANSFBUILDER_HPP_

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceHNDoFMapBase.hpp"
#include "SundanceMixedDOFMap.hpp"
#include "SundanceTransformationHN.hpp"
#include "SundanceTransformationBase.hpp"
#include "SundanceInequivalentElementTransformation.hpp"

#include "PlayaVectorType.hpp"
#include "PlayaVectorDecl.hpp"

namespace Sundance {

  using namespace Teuchos;
  using namespace Playa;

  /** This builds a transformation for one discrete space <br>
   *  The class also calls the transformation method of the created transformation object
   *  */
  class DiscreteSpaceTransfBuilder {
  public:

    /** Empty Ctor, in this case there will be no transformation*/
    DiscreteSpaceTransfBuilder();

    /** */
    DiscreteSpaceTransfBuilder( const Mesh& mesh, const BasisArray& basis,
				const RCP<DOFMapBase>& map );

    /** Dtor */
    virtual ~DiscreteSpaceTransfBuilder() {;}

    /** If there is a valid transformation then returns true, else false */
    const inline bool validTransformation() const { return hasTransformation_; }

	/** Function to make the local transformation for the discrete space
	 * @param dofs                 [in] the array with the DoF numbers
	 * @param functionIDs           [in] the functionIDs of this chunk
	 * @param chunkID              [in] function chunk ID
	 * @param nrDoFsPerCell        [in] nr DoF per Cell (per one chunk)
	 * @param nrFunctions          [in] nr of functions in this chunk
	 * @param ghostView            [in] the array where we have to get values from
	 * @param localValues          [out]  the transformed values */
	void getDoFsWithTransformation(const Array<int>& dofs,
			                             const Array<int>& functionIDs ,
			                             const int chunkID ,
			                             const int nrDoFsPerCell ,
			                             const int nrFunctions ,
			                 			 const int cellDim ,
			                 			 const Array<int>& cellLIDs,
			                             const RCP<GhostView<double> >& ghostView ,
			                             Array<double>& localValues) const;

  protected:

    /** for verbosity */
    int verb() const {return verb_;}

  private:

    /** verbosity */
    int verb_;

    /** Number of different function defined on the DoFMap*/
    int basisSize_;

    /** true if has a valid transformation */
    mutable bool hasTransformation_;

    /** The transformation for the Discrete space */
    mutable RCP<TransformationBase> transformation_;

  };

}

#endif /* SUNDANCEDISCRETESPACETRANSFBUILDER_HPP_ */
