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

#ifndef SUNDANCE_FIXEDARRAY_H
#define SUNDANCE_FIXEDARRAY_H


#include "SundanceDefs.hpp"
#include <iostream>
#include <algorithm>

namespace Sundance
{

/**
 * A simple fixed-size array 
 */
template <int N, class T> class FixedArray
{
public:
  /** */
  FixedArray(){;}

  /** */
  const T& operator[](int i) const 
    {
      TEUCHOS_TEST_FOR_EXCEPT(i<0);
      TEUCHOS_TEST_FOR_EXCEPT(i>=N);
      return data_[i];
    }

  /** */
  T& operator[](int i) 
    {
      TEUCHOS_TEST_FOR_EXCEPT(i<0);
      TEUCHOS_TEST_FOR_EXCEPT(i>=N);
      return data_[i];
    }

  /** */
  int size() const {return N;}

  /** */
  const T* begin() const {return &(data_[0]);}

  /** */
  const T* end() const {return begin()+N;}

  /** */
  T* begin() {return data_;}
  
  /** */
  T* end() {return data_+N;}

  /** */
  bool operator<(const FixedArray<N, T>& other) const
    {
      return std::lexicographical_compare(begin(), end(), 
        other.begin(), other.end());
    }

private:
  T data_[N];
};



/** \relate FixedArray */
template<class T>
FixedArray<2,T> makeFA(const T& a, const T& b)
{
  FixedArray<2,T> rtn;
  rtn[0] = a;
  rtn[1] = b;
  return rtn;
}

/** \relate FixedArray */
template<class T>
FixedArray<3,T> makeFA(const T& a, const T& b, const T& c)
{
  FixedArray<3,T> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  return rtn;
}


/** \relate FixedArray */
template<class T>
FixedArray<4,T> makeFA(const T& a, const T& b, const T& c, const T& d)
{
  FixedArray<4,T> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[2] = d;
  return rtn;
}

}


namespace std
{
/** \relates FixedArray */
template <int N, class T> 
ostream& operator<<(std::ostream& os, const Sundance::FixedArray<N, T>& f)
{
  os << "{";
  for (int i=0; i<N; i++) 
  {
    if (i>0) os << ", ";
    os << f[i];
  }
  os << "}";
  return os;
}

}

#endif
