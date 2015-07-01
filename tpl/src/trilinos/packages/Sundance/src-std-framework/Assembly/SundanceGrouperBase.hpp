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

#ifndef SUNDANCE_GROUPERBASE_H
#define SUNDANCE_GROUPERBASE_H

#include "SundanceDefs.hpp"
#include "SundanceCellType.hpp"
#include "SundanceParametrizedCurve.hpp"
#include "SundanceMesh.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Sundance
{
using namespace Teuchos;

class EquationSet;
class SparsitySuperset;
class MultiIndex;
class MultipleDeriv;
class QuadratureFamily;
class BasisFamily;
class IntegralGroup;


/** 
 * Grouper
 */
class GrouperBase
{
public:
  /** */
  GrouperBase() {}

  /** */
  virtual ~GrouperBase(){;}

  /** */
  virtual void findGroups(const EquationSet& eqn,
    const CellType& maxCellType,
    int spatialDim,
    const CellType& cellType,
    int cellDim,
    const QuadratureFamily& quad,
    const RCP<SparsitySuperset>& sparsity,
    bool isInternalBdry,
    Array<RCP<IntegralGroup> >& groups,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh ) const = 0 ;


  /** */
  void setVerb(
    int setupVerb,
    int integrationVerb,
    int transformVerb);
    

  /** */
  int setupVerb() const {return setupVerb_;}
    
  /** */
  int integrationVerb() const {return integrationVerb_;}
    
  /** */
  int transformVerb() const {return transformVerb_;}

protected:
  void extractWeakForm(const EquationSet& eqn,
    const MultipleDeriv& functionalDeriv,
    BasisFamily& testBasis, 
    BasisFamily& unkBasis,
    MultiIndex& miTest, MultiIndex& miUnk,
    int& rawVarID, int& rawUnkID,  
    int& reducedTestID, int& reducedUnkID, 
    int& testBlock, int& unkBlock, 
    int& rawParamID, int& reducedParamID,
    bool& isOneForm, bool& hasParam) const ;
                              
private:
  int setupVerb_;
  int integrationVerb_;
  int transformVerb_;
};

}


#endif
