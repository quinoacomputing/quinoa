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

#ifndef SUNDANCE_MULTIINDEX_H
#define SUNDANCE_MULTIINDEX_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_XMLObject.hpp"
#include <string>
#include <stdexcept>

namespace Sundance
{
using namespace Teuchos;

/**
 * An integer vector representing a multivariate derivative.
 */

class MultiIndex
{
public:
  /** constructs D(0,0,0) */
  MultiIndex();
  /** constructs a multiindex D(x,y,z) */
  MultiIndex(int x, int y, int z);

  /** */
  bool operator==(const MultiIndex& other) const ;

  /** */
  bool operator<(const MultiIndex& other) const ;

  /** */
  const int& operator[](int i) const {return m_[i];}

  /** */
  int& operator[](int i) {return m_[i];}

  /** */
  MultiIndex operator+(const MultiIndex& other) const ;

  /** */
  MultiIndex operator-(const MultiIndex& other) const ;

  /** */
  MultiIndex operator-() const ;

  /** */
  std::string toString() const ;

  /** */
  XMLObject toXML() const ;

  /** */
  int order() const ;

  /** */
  int firstOrderDirection() const ;

  /** */
  static int maxDim() {return 3;}

  /** */
  bool isValid() const ;

  /** */
  std::string coordForm() const ;
private:
  Array<int> m_;
};
}

namespace Teuchos
{

/** \relates Sundance::MultiIndex */
inline std::string toString(const Sundance::MultiIndex& h)
{return h.toString();}

}

namespace std
{
/** \relates Sundance::MultiIndex */
inline ostream& operator<<(std::ostream& os, 
  const Sundance::MultiIndex& h)
{
  os << h.toString();
  return os;
}
}

#endif
