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

#include "SundanceUserDefFunctor.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvalVector.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


UserDefFunctor::UserDefFunctor(const std::string& name, 
                               int domainDim, 
                               int rangeDim)
  :  name_(name), 
     elemNames_(), 
     domainDim_(domainDim), 
     rangeDim_(rangeDim)
{
  TEUCHOS_TEST_FOR_EXCEPT(domainDim_ <= 0);
  TEUCHOS_TEST_FOR_EXCEPT(rangeDim_ <= 0);

  elemNames_.resize(rangeDim_);

  if (rangeDim_==1) 
    {
      elemNames_[0] = name_;
    }
  else
    {
      for (int i=0; i<rangeDim_; i++)
        {
          elemNames_[i] = name_ + "[" + Teuchos::toString(i) + "]";
        }
    }
}





void UserDefFunctor::eval0(const Array<double>& in, double* out) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                     "eval0() function not supported for functor " << name_);
}

void UserDefFunctor::eval1(const Array<double>& in, double* out, double* outDerivs) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                     "eval1() function not supported for functor " << name_);
}

