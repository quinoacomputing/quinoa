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

#include "SundanceSparsitySuperset.hpp"
#include "SundanceEvaluatableExpr.hpp"

#include "SundanceEvalVector.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"



using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;




SparsitySuperset::SparsitySuperset(const Set<MultipleDeriv>& C,
                                   const Set<MultipleDeriv>& V)
  : maxOrder_(0),
    derivToIndexMap_(),
    derivs_(),
    states_(),
    multiIndex_(),
    numConstantDerivs_(0),
    numVectorDerivs_(0)
{
  for (Set<MultipleDeriv>::const_iterator i=C.begin(); i!=C.end(); i++)
    {
      addDeriv(*i, ConstantDeriv);
    }

  for (Set<MultipleDeriv>::const_iterator i=V.begin(); i!=V.end(); i++)
    {
      addDeriv(*i, VectorDeriv);
    }
}



void SparsitySuperset::addDeriv(const MultipleDeriv& d, 
                               const DerivState& state)
{
  maxOrder_ = std::max(d.order(), maxOrder_);

  if (containsDeriv(d))
    {
      const DerivState& oldState = states_[getIndex(d)];
      if (state > oldState) 
        {
          states_[getIndex(d)] = state;
          numConstantDerivs_--;
          numVectorDerivs_++;
        }
    }
  else
    {
      int index = derivs_.size();
      derivs_.append(d);
      states_.append(state);
      derivToIndexMap_.put(d, index);
      MultiIndex mi;
      for (MultipleDeriv::const_iterator i=d.begin(); 
           i != d.end(); i++)
        {
          if (i->isCoordDeriv())
            {
              MultiIndex m;
              int dir = i->coordDerivDir();
              m[dir] = 1;
              mi = mi + m;
            }

        }
      multiIndex_.append(mi);
      if (state==ConstantDeriv) numConstantDerivs_++;
      else numVectorDerivs_++;
    }
}

void SparsitySuperset::addDeriv(const Deriv& d, 
                               const DerivState& state)
{
  MultipleDeriv md;
  md.put(d);
  addDeriv(md, state);
}



bool SparsitySuperset::containsDeriv(const MultipleDeriv& d) const
{
  return derivToIndexMap_.containsKey(d);
}

int SparsitySuperset::getIndex(const MultipleDeriv& d) const
{
  if (!containsDeriv(d)) return -1;
  return derivToIndexMap_.get(d);
}



void SparsitySuperset::print(std::ostream& os,
                             const Array<RCP<EvalVector> >& vecResults,
                             const Array<double>& constantResults) const 
{
  Tabs tabs;

  /* find the maximum size of the std::string reps of the derivatives.
   * We'll use this to set the field width for printing derivatives. */
  int maxlen = 25;
  for (int i=0; i<derivs_.size(); i++)
    {
      int s = derivs_[i].toString().length();
      if (s > maxlen) maxlen = s;
    }


  int vecIndex=0;
  int constIndex = 0;
  os << tabs << "Results Superset" << std::endl;
  for (int i=0; i<derivs_.size(); i++)
    {
      os << tabs << i << "\t\t";
      os.width(maxlen);
      os.setf(std::ios_base::left, std::ios_base::adjustfield);
      os << derivs_[i].toString() << "\t\t" ;
      switch(states_[i])
        {
        case ZeroDeriv:
          os  << "Zero" << std::endl;
          break;
        case ConstantDeriv:
          os << constantResults[constIndex++] << std::endl;
          break;
        case VectorDeriv:
          if (vecResults[vecIndex].get()==0)
            {
              os << "{Null}";
            }
          else
            {
              vecResults[vecIndex]->print(os);
            }
          vecIndex++;
          os << std::endl;
          break;
        }
    }
}



void SparsitySuperset::print(std::ostream& os) const 
{
  Tabs tabs;

  /* find the maximum size of the std::string reps of the derivatives.
   * We'll use this to set the field width for printing derivatives. */
  int maxlen = 25;
  for (int i=0; i<derivs_.size(); i++)
    {
      int s = derivs_[i].toString().length();
      if (s > maxlen) maxlen = s;
    }


  os << tabs << "SparsitySuperset" << std::endl;
  for (int i=0; i<derivs_.size(); i++)
    {
      os << tabs << i << "\tderiv=\t";
      os.width(maxlen);
      os.setf(std::ios_base::left, std::ios_base::adjustfield);
      os << derivs_[i].toString() << "\tstate=\t" ;
      switch(states_[i])
        {
        case ZeroDeriv:
          os  << "Zero" << std::endl;
          break;
        case ConstantDeriv:
          os << "Constant" << std::endl;
          break;
        case VectorDeriv:
          os << "Vector" << std::endl;
          break;
        }
    }
}

string SparsitySuperset::toString() const 
{
	std::ostringstream ss;
	print(ss);
	return ss.str();
}

DerivSet SparsitySuperset::derivSet() const
{
  DerivSet rtn;
  for (int i=0; i<numDerivs(); i++) 
    {
      if (state(i) != ZeroDeriv) rtn.put(deriv(i));
    }
  return rtn;
}
