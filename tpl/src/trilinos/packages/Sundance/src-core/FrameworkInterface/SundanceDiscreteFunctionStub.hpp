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

#ifndef SUNDANCE_DISCRETEFUNCTIONSTUB_H
#define SUNDANCE_DISCRETEFUNCTIONSTUB_H

#include "SundanceDefs.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceDiscreteFuncDataStub.hpp"
#include "SundanceSpectralBasis.hpp"

namespace Sundance
{

/** 
 * DiscreteFunctionStub is the base class for discrete functions. 
 * Each framework will need to implement its own subclass of
 * DiscreteFunctionStub. 
 *
 * The interface is left very minimal so as to not place
 * any constraints on how a framework might specify vectors
 * and bases. When a framework needs any information about the
 * discrete function, it will have to get it by downcasting
 * to the appropriate framework-specific subclass.
 *
 * <h4> Writing a DiscreteFunctionStub subclass </h4>
 *
 * For purposes of interaction with the Sundance core, no 
 * additional methods are required.
 * However, most frameworks will require extensions to 
 * DiscreteFunctionStub that can supply the framework with information
 * on the basis and vector used by the discrete func. See the
 * demo and standard frameworks for information on how to do this.
 */
class DiscreteFunctionStub : public ListExpr
{
public:
  /** */
  DiscreteFunctionStub(const std::string& name, 
    int tensorOrder=0, int dim=1,
    const RCP<DiscreteFuncDataStub>& data=RCP<DiscreteFuncDataStub>(),
    int listIndex=0);
  /** */
  DiscreteFunctionStub(const Array<string>& name, 
    const Array<std::pair<int,int> >& tensorStructure,
    const RCP<DiscreteFuncDataStub>& data=RCP<DiscreteFuncDataStub>());
     
  /** */
  DiscreteFunctionStub(const std::string& name, 
    const SpectralBasis& sbasis, 
    int tensorOrder=0, int dim=1,
    const RCP<DiscreteFuncDataStub>& data=RCP<DiscreteFuncDataStub>(),
    int listIndex=0);
     
  /** */
  DiscreteFunctionStub(const Array<string>& name, 
    const SpectralBasis& sbasis, 
    const Array<std::pair<int,int> >& tensorStructure,
    const RCP<DiscreteFuncDataStub>& data=RCP<DiscreteFuncDataStub>());

  /** virtual destructor */
  virtual ~DiscreteFunctionStub() {;}

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

protected:
  /** */
  void initTensor(const std::string& name, 
    int tensorOrder, int dim,
    const RCP<DiscreteFuncDataStub>& data,
    int listIndex);
  /** */
  void initTensorSpectral(const std::string& name, 
    const SpectralBasis& sbasis, 
    int tensorOrder, int dim,
    const RCP<DiscreteFuncDataStub>& data,
    int listIndex);

  /** */
  const RCP<DiscreteFuncDataStub>& dataStub() const {return data_;}

private:
  RCP<DiscreteFuncDataStub> data_;
};
}



#endif
