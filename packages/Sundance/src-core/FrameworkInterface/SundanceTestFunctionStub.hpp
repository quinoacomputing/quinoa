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

#ifndef SUNDANCE_TESTFUNCTIONSTUB_H
#define SUNDANCE_TESTFUNCTIONSTUB_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceTestFuncDataStub.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"

namespace Sundance
{
class DiscreteFunctionStub;

using namespace Teuchos;

/** 
 * TestFunctionStub is the base class for unknown functions. 
 * Each framework will need to implement its own subclass of
 * TestFunctionStub. 
 *
 * The interface is left very minimal so as to not place
 * any constraints on how a framework might specify the basis.
 * When a framework needs any information about the
 * unknown function, it will have to get it by downcasting
 * to the appropriate framework-specific subclass.
 *
 * <h4> Writing a TestFunctionStub subclass </h4>
 *
 * For purposes of interaction with the Sundance core, no 
 * additional methods are required.
 * However, most frameworks will require extensions to 
 * TestFunctionStub that can supply the framework with information
 * on the basis used by the unknown func. See the
 * demo and standard frameworks for information on how to do this.
 */
class TestFunctionStub : public SymbolicFunc
{
public:
  /** */
  TestFunctionStub(const std::string& name, 
    int tensorOrder=0, int dim=1,
    const RCP<const TestFuncDataStub>& data=RCP<const TestFuncDataStub>());

  /** */
  TestFunctionStub(const std::string& name, const SpectralBasis& sbasis, 
    int tensorOrder=0, int dim=1,
    const RCP<const TestFuncDataStub>& data=RCP<const TestFuncDataStub>());

  /** virtual destructor */
  virtual ~TestFunctionStub() {;}

  /** */
  bool isTestFunction() const {return true;}

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

protected:

  /** */
  const RCP<const TestFuncDataStub>& dataStub() const {return data_;}

private:
  RCP<const TestFuncDataStub> data_;


};

}
                  

#endif
