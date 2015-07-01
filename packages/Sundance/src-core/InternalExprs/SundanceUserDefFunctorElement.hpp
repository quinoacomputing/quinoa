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

#ifndef SUNDANCE_USERDEFFUNCTORELEMENT_H
#define SUNDANCE_USERDEFFUNCTORELEMENT_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceMultiSet.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "Teuchos_Array.hpp"




namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




/**
 * Scalar-valued element of a vector-valued functor
 */
class UserDefFunctorElement
{
public:
  /** ctor */
  UserDefFunctorElement(const RCP<const UserDefFunctor>& functor,
    int myIndex);

  /** */
  virtual ~UserDefFunctorElement(){;}

  /** */
  const std::string& name() const {return master_->name(myIndex());}

  /** */
  const std::string& masterName() const {return master_->name();}

  /** */
  void evalArgDerivs(int maxOrder, 
    const Array<double>& in,
    Array<double>& outDerivs) const ;

  /** */
  void getArgDerivIndices(const Array<int>& orders,
    Sundance::Map<MultiSet<int>, int>& varArgDerivs,
    Sundance::Map<MultiSet<int>, int>& constArgDerivs) const ;

  /** */
  int numArgs() const {return master_->domainDim();}

  /** */
  void reset() const {master_->reset();}

  /** Return the index of this element into the list-valued 
   * user defined op */
  int myIndex() const {return myIndex_;}

  /** */
  const UserDefFunctor* master() const 
    {return master_.get();}

  /** */
  int maxOrder() const {return master_->maxOrder();}


private:
  const RCP<const UserDefFunctor> master_;
  const int myIndex_;
};


}


#endif
