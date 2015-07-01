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

#ifndef SUNDANCE_ALGEBRASPECIFIER_HPP
#define SUNDANCE_ALGEBRASPECIFIER_HPP

#include "SundanceEnumTypeField.hpp"

namespace Sundance
{
using namespace Sundance;

/** */
enum AlgebraType {VectorAT, NormalAT, CoordCompAT, ScalarAT};

/** AlgebraSpecifier is used to indicate whether a given object represents
 * a vector, scalar, component of a vector, etc. 
 */
class AlgebraSpecifier : public EnumTypeField<AlgebraType>
{
public:
  /** */
  AlgebraSpecifier();

  /** */
  AlgebraSpecifier(int direction);

  /** */
  AlgebraSpecifier(const AlgebraType& vct);

  /** Return true iff I am a scalar */
  bool isScalar() const {return isType(ScalarAT);}

  /** Return true iff I am a normal component of a vector */
  bool isNormal() const {return isType(NormalAT);}

  /** Return true iff I am a vector */
  bool isVector() const {return isType(VectorAT);}

  /** Return true iff I am a coordinate component of a vector */
  bool isCoordinateComponent() const {return isType(CoordCompAT);}

  /** Return coordinate direction if I am a coordinate component, 
   * otherwise, throw an exception. */
  int direction() const ;

  /** Write a std::string representation */
  std::string toString() const ;

  /** Comparison operator for use in sorted containers */
  bool operator<(const AlgebraSpecifier& other) const ;

private:
  int direction_;
};

/** \relates AlgebraSpecifier */
AlgebraSpecifier vectorAlgebraSpec();

/** \relates AlgebraSpecifier */
AlgebraSpecifier scalarAlgebraSpec();

/** \relates AlgebraSpecifier */
AlgebraSpecifier normalAlgebraSpec();

/** \relates AlgebraSpecifier */
AlgebraSpecifier coordAlgebraSpec(int dir);

}


namespace std
{
/** \relates AlgebraSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::AlgebraSpecifier& vcs);
/** \relates AlgebraSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::AlgebraType& vct);
}


#endif
