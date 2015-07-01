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

#ifndef SUNDANCE_TESTFUNCELEMENT_H
#define SUNDANCE_TESTFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceTestFuncDataStub.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




/** 
 * TestFuncElement represents a scalar-valued element of a (possibly)
 * list-valued TestFunction
 */
class TestFuncElement : public SymbolicFuncElement
{
public:
  /** */
  TestFuncElement(const RCP<const TestFuncDataStub>& commonData,
    const std::string& name,
    const std::string& suffix, const FunctionIdentifier& fid);

  /** virtual destructor */
  virtual ~TestFuncElement() {;}

  /** Test whether all terms have test functions. 
   * I'm a test function, so return true */
  virtual bool everyTermHasTestFunctions() const {return true;}

  /** Test whether this expr contains a test function. 
   * I'm a test function, so return true. */
  virtual bool hasTestFunctions() const {return true;}

  /** */
  virtual bool isTestFunction() const {return true;}

  /** 
   * Indicate whether the expression is linear 
   * with respect to test functions */
  virtual bool isLinearInTests() const {return true;}
      
  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;


  /** */
  virtual XMLObject toXML() const ;

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}
      
private:
};
}

#endif
