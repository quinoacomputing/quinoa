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

#ifndef SUNDANCE_SPATIALDERIVSPECIFIER_HPP
#define SUNDANCE_SPATIALDERIVSPECIFIER_HPP

#include "SundanceEnumTypeField.hpp"
#include "SundanceMultiIndex.hpp"

namespace Sundance
{
using namespace Sundance;

/** */
enum SpatialDerivType {IdentitySDT, PartialSDT, NormalSDT, DivSDT};



/** 
 * This class is a compact description of type of spatial derivative
 * acting on an operative function: partial derivative, divergence, 
 * or normal derivative. 
 */
class SpatialDerivSpecifier : public EnumTypeField<SpatialDerivType>
{
public:
  /** Empty ctor creates an identity operator 
   * (zeroth order partial derivative) */
  SpatialDerivSpecifier();

  /** Create a spatial derivative */
  SpatialDerivSpecifier(const MultiIndex& mi);

  /** Create a derivative of a specified type and order. */
  SpatialDerivSpecifier(const SpatialDerivType& type, int order=0);

  /** Return the multiindex of a spatial partial derivative */
  const MultiIndex& mi() const ;

  /** Return true if I am a divergence */
  bool isDivergence() const ;

  /** Return true if I am a partial derivative in a coordinate direction */
  bool isPartial() const ;

  /** Return true if I am a normal derivative */
  bool isNormal() const ;

  /** Return true if I am an identity operator */
  bool isIdentity() const ;

  /** Return the order of differentiation in the normal direction */
  int normalDerivOrder() const ;

  /** Return the order of differentiation */
  int derivOrder() const ;

  /** Write me to a std::string */
  std::string toString() const ;

  /** Comparison operator for use in sorted containers */
  bool operator<(const SpatialDerivSpecifier& other) const ;

  /** Create a new derivative that increments my multiindex by the input
   * multiindex */
  SpatialDerivSpecifier derivWrtMultiIndex(const MultiIndex& mi) const ;

private:
  MultiIndex mi_;

  int normalDerivOrder_;
};


}


namespace std
{
/** \relates SpatialDerivSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::SpatialDerivSpecifier& sds);
/** \relates SpatialDerivSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::SpatialDerivType& sdt);
}


#endif
