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

#ifndef SUNDANCE_FUNCTIONALGRADIENTASSEMBLYKERNEL_H
#define SUNDANCE_FUNCTIONALGRADIENTASSEMBLYKERNEL_H

#include "SundanceDefs.hpp"
#include "SundanceVectorFillingAssemblyKernel.hpp"
#include "SundanceFunctionalAssemblyKernel.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * FunctionalGradientAssemblyKernel does assembly of a functional and
 * its gradient. 
 */
class FunctionalGradientAssemblyKernel : public AssemblyKernelBase
{
public:
  /** */
  FunctionalGradientAssemblyKernel(const MPIComm& comm,
    const Array<RCP<DOFMapBase> >& dofMap,
    const Array<RCP<Array<int> > >& isBCIndex,
    const Array<int>& lowestLocalIndex,
    Array<Vector<double> >& grad,
    bool partitionBCs,
    double* value, 
    int verb)
    : AssemblyKernelBase(verb),
      funcKernel_(rcp(new FunctionalAssemblyKernel(comm, value, verb))),
      vecKernel_(rcp(new VectorAssemblyKernel(dofMap, isBCIndex,
            lowestLocalIndex, grad, partitionBCs, verb)))
    {
      TEUCHOS_TEST_FOR_EXCEPTION(grad.size() != 1, std::logic_error,
        "assembly target in FunctionalGradientAssemblyKernel should not "
        "be a multivector");
    }

  /** */
  void prepareForWorkSet(
    const Array<Set<int> >& requiredTests,
    const Array<Set<int> >& requiredUnks,
    RCP<StdFwkEvalMediator> mediator) 
    {
      funcKernel_->prepareForWorkSet(requiredTests, requiredUnks, mediator);
      vecKernel_->prepareForWorkSet(requiredTests, requiredUnks, mediator);
    }

  /** */
  void fill(bool isBC,
    const IntegralGroup& group,
    const RCP<Array<double> >& localValues) 
    {
      if (group.isOneForm())
      {
        vecKernel_->fill(isBC, group, localValues);
      }
      else if (group.isZeroForm())
      {
        funcKernel_->fill(isBC, group, localValues);
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPT(group.isTwoForm());
      }
    }

  /** */
  void postLoopFinalization()
    {
      funcKernel_->postLoopFinalization();
      vecKernel_->postLoopFinalization();
    }

  /** */
  void setVerb(int verb) 
    {
      AssemblyKernelBase::setVerb(verb);
      funcKernel_->setVerb(verb);
      vecKernel_->setVerb(verb);
    }

private:
  RCP<FunctionalAssemblyKernel> funcKernel_;
  RCP<VectorAssemblyKernel> vecKernel_;
};


}



#endif
