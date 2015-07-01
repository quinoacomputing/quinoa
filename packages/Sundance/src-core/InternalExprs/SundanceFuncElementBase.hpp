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

#ifndef SUNDANCE_FUNCELEMENTBASE_H
#define SUNDANCE_FUNCELEMENTBASE_H


#include "SundanceDefs.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceFunctionWithID.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;





/** 
 * FuncElementBase defines the interface for scalar-valued elements
 * of Sundance functions. At the user level, Sundance functions can be
 * list (e.g, vector or tensor) valued; internally, however, compound
 * expressions use only scalar functions deriving from the 
 * FuncElementBase class. 
 */
class FuncElementBase : public virtual ScalarExpr,
  public FunctionWithID
{
public:
  /** */
  FuncElementBase(const std::string& rootName,
    const std::string& suffix,
    const FunctionIdentifier& fid);
  /** */
  FuncElementBase(const std::string& rootName);

  /** virtual destructor */
  virtual ~FuncElementBase() {;}

  /** Return the name of this function */
  const std::string& name() const {return name_;}

  /** Return the root name of this function */
  const std::string& rootName() const {return rootName_;}

  /** Return the root name of this function */
  const std::string& suffix() const {return suffix_;}

  /** Write self in text form */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;


protected:
private:

  std::string name_;

  std::string rootName_;

  std::string suffix_;
};

}

#endif
