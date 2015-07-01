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

#ifndef SUNDANCE_EXPRBASE_H
#define SUNDANCE_EXPRBASE_H


#include "SundanceDefs.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include "PlayaHandleable.hpp"



namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;





/** */
class ExprBase : public Playa::Handleable<ExprBase>
{
public:
  /** empty ctor */
  ExprBase();

  /** virtual destructor */
  virtual ~ExprBase() {;}

  /** Write a simple text description suitable 
   * for output to a terminal */
  virtual std::ostream& toText(std::ostream& os, bool paren) const = 0 ;

  /** Append to the set of func IDs present in this expression.
   * Base class does nothing */
  virtual void accumulateFuncSet(Set<int>& funcIDs, 
    const Set<int>& activeSet) const {;}

  /** Indicate whether this expression contains any test 
   * functions. Default is to return false. This will be
   * overridden by TestFuncElement and ExprWithChildren. */
  virtual bool hasTestFunctions() const {return false;}
  /** 
   * Indicate whether the expression contains unknown functions */
  virtual bool hasUnkFunctions() const {return false;}

  /** */
  std::string toString() const ;

  /** Write in XML */
  virtual XMLObject toXML() const = 0 ;

  /** Return a descriptive name for the expression subtype */
  virtual std::string typeName() const ;

protected:
};



}

#endif
