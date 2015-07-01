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
 * SundanceAssemblyTransformationBuilder.hpp
 *
 *  Created on: Mar 16, 2010
 *      Author: benk
 */

#ifndef SUNDANCEASSEMBLYTRANSFORMATIONBUILDER_HPP_
#define SUNDANCEASSEMBLYTRANSFORMATIONBUILDER_HPP_

#include "SundanceIntegralGroup.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceCellType.hpp"
#include "SundanceTransformationBase.hpp"
#include "SundanceTransformationHN.hpp"
#include "SundanceInequivalentElementTransformation.hpp"

namespace Sundance {

  class AssemblyTransformationBuilder {

  public:


    /** This is (mainly used from the Assemby Loop) <br>
	In case when the basis and the test space are different <br>
	An Integral group migh have many members, but the test and trial space combination is UNIQUE !!!
	Hanging node mesh with Hermit base might be tricky ... */
    AssemblyTransformationBuilder( const RCP<IntegralGroup>& group ,
				   const Array<RCP<DOFMapBase> >& rowMaps ,
				   const Array<RCP<DOFMapBase> >& colMaps ,
				   const Mesh& mesh);

    /** */
    virtual ~AssemblyTransformationBuilder();


    /** */
    void applyTransformsToAssembly( int groupIndex ,
				    int entryPerCell ,
				    CellType cellType ,
				    int cellDim,
				    CellType maxCellType ,
				    const CellJacobianBatch& JTrans,
				    const CellJacobianBatch& JVol,
				    const Array<int>& facetNum,
				    const RCP<Array<int> >& cellLIDs,
				    RCP<Array<double> >& A );


    /** return the preTransformation*/
    const RCP<TransformationBase>& getPreTransformation() const { return preTransformation_; }

    /** return the postTransformation*/
    const RCP<TransformationBase>& getPostTransformation() const { return postTransformation_; }

    /** verbosity level */
    int verb() const {return verb_;}

    /** set verbosity level */
    void setVerb(int verb) {verb_=verb;}

  private:

    /* verbosity */
    int verb_;

    /* number of columns of the matrix/vector which should be transformed */
    int nrCol_;

    /* number of rows of the matrix/vector which should be transformed */
    int nrRow_;

    /** The pre-transformation */
    mutable RCP<TransformationBase> preTransformation_;

    /** The post-transformation (might be the same as the pre-transformation, the RCP will take care of that) */
    mutable RCP<TransformationBase> postTransformation_;

    /** */
    int testFuncID_;

    /** */
    int unkFuncID_;

    /** */
    const DOFMapBase* _myRowDOFMap;

    /** */
    const DOFMapBase* _myColDOFMap;

    /** */
    mutable bool hasTransformation_;

    /** */
    bool hasPreTransformation_;
    bool hasPostTransformation_;

    /** */
    mutable bool onlyVectorTransformation_;
  };

}

#endif /* SUNDANCEASSEMBLYTRANSFORMATIONBUILDER_HPP_ */
