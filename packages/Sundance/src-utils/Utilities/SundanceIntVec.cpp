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


#include "SundanceIntVec.hpp"
#include "SundanceDebug.hpp"
#include "SundanceCombinatorialUtils.hpp"

namespace Sundance
{
using Teuchos::Array;

IntVec::IntVec(int n)
  : data_(n)
{
  for (int i=0; i<n; i++) data_[i] = 0;
}

IntVec IntVec::operator+(const IntVec& other) const
{
  TEUCHOS_TEST_FOR_EXCEPT(size() != other.size());

  IntVec rtn(size());
  for (int i=0; i<size(); i++) rtn[i] = data_[i] + other[i];
  
  return rtn;
}


IntVec IntVec::operator*(int a) const
{
  IntVec rtn(size());
  for (int i=0; i<size(); i++) rtn[i] = a*data_[i];
  
  return rtn;
}

int IntVec::factorial() const
{
  int rtn=1;

  for (int i=0; i<size(); i++)
  {
    int n_i = data_[i];
    for (int j=1; j<=n_i; j++) rtn *= j;
  }
  return rtn;
}


int IntVec::pow(const IntVec& other) const
{
  TEUCHOS_TEST_FOR_EXCEPT(size() != other.size());
  int rtn=1;

  for (int i=0; i<size(); i++)
  {
    int n_i = data_[i];
    int p_i = other[i];
    for (int j=1; j<=p_i; j++) rtn *= n_i;
  }
  return rtn;
}

int IntVec::abs() const
{
  int rtn = 0;
  for (int i=0; i<size(); i++)
  {
    rtn += ::abs(data_[i]);
  }
  return rtn;
}

int IntVec::norm() const
{
  int rtn = 0;
  for (int i=0; i<size(); i++)
  {
    if (::abs(data_[i]) > rtn) rtn = ::abs(data_[i]);
  }
  return rtn;
}

bool IntVec::operator==(const IntVec& other) const
{
  if (size() != other.size()) return false;
  for (int i=0; i<size(); i++)
  {
    if (data_[i] != other.data_[i]) return false;
  }
  return true;
}

bool IntVec::operator<(const IntVec& other) const
{
  if (size() < other.size()) return true;
  if (size() > other.size()) return false;
  if (abs() < other.abs()) return true;
  if (abs() > other.abs()) return false;

  for (int i=0; i<size(); i++)
  {
    if (data_[i] < other.data_[i]) return true;
    if (data_[i] > other.data_[i]) return false;
  }

  return false;
}

void IntVec::print(std::ostream& os) const 
{
  os << "IVec[" << data_ << "]";
}


void IntVec::getPartitions(int M, Array<Array<IntVec> >& parts) const
{
  Array<Array<Array<int> > > rComp(size());
  Array<int> radix(size());
  
  for (int i=0; i<size(); i++)
  {
    restrictedCompositions(data_[i], M, rComp[i]);
    radix[i] = rComp[i].size();
  }

  Array<int> pick(size(), 0);
  bool workLeft = true;
  while (workLeft)
  {
    Array<IntVec> part(M);
    for (int j=0; j<M; j++) part[j] = IntVec(size());
    for (int i=0; i<size(); i++)
    {
      TEUCHOS_TEST_FOR_EXCEPT(pick[i] >= rComp[i].size());
      const Array<int>& p = rComp[i][pick[i]];
      TEUCHOS_TEST_FOR_EXCEPT(p.size() != M);
      for (int j=0; j<M; j++) part[j][i] = p[j];
    }
    bool isNonzero = true;
    for (int i=0; i<part.size(); i++) 
    {
      if (part[i].abs()==0) {isNonzero = false; break;}
    }
    if (isNonzero) parts.append(part);
    workLeft = nextNum(pick, radix);
  }
}

}



	





