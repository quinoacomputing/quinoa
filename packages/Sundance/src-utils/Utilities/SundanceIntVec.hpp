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

#ifndef SUNDANCE_INTVEC_H
#define SUNDANCE_INTVEC_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "PlayaPrintable.hpp"

namespace Sundance
{
using Teuchos::Array;
using Teuchos::tuple;

/** 
 * An integer vector class for use in the generalized chain rule
 * of Constantine and Savits (1996). See the paper by CS for 
 * definitions of the various operations.
 */
class IntVec : public Playa::Printable
{
public:
  /** */
  IntVec(){}
  /** */
  IntVec(const Array<int>& d) : data_(d) {}
  /** */
  IntVec(int n);

  /** */
  int size() const {return data_.size();}

  /** */
  int operator[](int i) const {return data_[i];}

  /** */
  int& operator[](int i) {return data_[i];}

  /** */
  IntVec operator+(const IntVec& other) const ;

  /** */
  IntVec operator*(int alpha) const ;

  /** Return the factorial of this vector, as defined by CS */
  int factorial() const ;

  /** */
  int pow(const IntVec& other) const ;

  /** Return the sum of elements in this vector */
  int abs() const ;

  /** Return the infinity norm of this vector */
  int norm() const ;

  /** */
  bool operator==(const IntVec& other) const ;

  /** */
  bool operator<(const IntVec& other) const ;

  /** */
  void print(std::ostream& os) const ;

  /** Get the length-M partitions of this vector. These are all 
   * list of exactly M vectors \f$ v_i\f$ such that
   * \f[ \sum_{i=1}^M v_i = this. \f]
   */
  void getPartitions(int M, Array<Array<IntVec> >& parts) const ;

private:
  Array<int> data_;
};



/** \relates IntVec */
inline IntVec operator*(int a, const IntVec& x)
{
  return x * a;
}

/** \relates IntVec*/
inline std::ostream& operator<<(std::ostream& os, const IntVec& v)
{
  v.print(os);
  return os;
}

/** */
inline IntVec intVec(int a)
{
  Array<int> dat = tuple(a);
  return dat;
}

/** */
inline IntVec intVec(int a, int b)
{
  Array<int> dat = tuple(a, b);
  return dat;
}

/** */
inline IntVec intVec(int a, int b, int c)
{
  Array<int> dat = tuple(a, b, c);
  return dat;
}

/** */
inline IntVec intVec(int a, int b, int c, int d)
{
  Array<int> dat = tuple(a, b, c, d);
  return dat;
}

/** */
inline IntVec intVec(int a, int b, int c, int d, int e)
{
  Array<int> dat = tuple(a, b, c, d, e);
  return dat;
}



}




#endif
