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


#include "SundanceEvalVectorArray.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;

EvalVectorArray::EvalVectorArray()
  : Array<RCP<EvalVector> >()
{;}

void EvalVectorArray::copy(const RCP<EvalVectorArray>& other)
{
  resize(other->size());

  for (int i=0; i<other->size(); i++)
    {
      (*this)[i]->copy((*other)[i]);
    }
}

void EvalVectorArray::steal(const RCP<EvalVectorArray>& other)
{
  resize(other->size());

  for (int i=0; i<other->size(); i++)
    {
      (*this)[i] = (*other)[i];
    }
}

ostream& EvalVectorArray::print(std::ostream& os, 
                                const SparsitySuperset* derivs) const
{
  Tabs tab;
  TEUCHOS_TEST_FOR_EXCEPTION(derivs->numDerivs() != this->size(),
                     std::logic_error,
                     "mismatch between deriv set size=" << derivs->numDerivs()
                     << "and result vector size " << this->size()
                     << "in EvalVectorArray::print");

  int maxlen = 25;
  for (int i=0; i<derivs->numDerivs(); i++)
    {
      int s = derivs->deriv(i).toString().length();
      if (s > maxlen) maxlen = s;
    }
  
  for (int i=0; i<derivs->numDerivs(); i++)
    {
      os << tab;
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << derivs->deriv(i).toString() << "\t\t";
      (*this)[i]->print(os);
      os << std::endl;
    }
  return os;
}
