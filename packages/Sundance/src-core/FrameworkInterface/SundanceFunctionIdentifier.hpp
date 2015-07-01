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

#ifndef SUNDANCE_FUNCTIONIDENTIFIER_H
#define SUNDANCE_FUNCTIONIDENTIFIER_H

#include "SundanceAlgebraSpecifier.hpp"

namespace Sundance
{

/** 
 * FunctionIdentifier provides a means for distinguishing between different
 * functions and different vector components of the same function. 
 * Functions discretized with vector bases will shared a common dofID,
 * because their vector components are not independent. Functions 
 * discretized componentwise will have different IDs for each component.
 */
class FunctionIdentifier
{
public:
  /** */
  FunctionIdentifier();
  /** ctor */
  FunctionIdentifier(const AlgebraSpecifier& algSpec);
  /** ctor */
  FunctionIdentifier(const FunctionIdentifier* parent,
    const AlgebraSpecifier& componentAlgSpec);

  /** */
  std::string toString() const ;

  /** Return the ID number to be used when assigning DOFs 
      for this function */
  int dofID() const {return dofID_;}

  /** If this FID corresponds to a vector component, return the 
   *  index of the coordinate direction */
  int componentIndex() const ;

  /** Return a specification of the type of object represented, i.e.,
   * a component in a coord direction, a normal component, or a whole
   * vector. */
  const AlgebraSpecifier& algSpec() const 
    {return algSpec_;}

  /** Create a new FID representing a component of "this" vector function. */
  FunctionIdentifier createComponent(int index) const ;

  /** Create a new FID representing the normal 
   * component of "this" vector function. */
  FunctionIdentifier createNormal() const ;

  /** Comparison operator for storage in sets and maps */
  bool operator<(const FunctionIdentifier& other) const ;

  /** Equality test */
  bool operator==(const FunctionIdentifier& other) const 
    {return !(*this!=other);} 

  /** Inequality test */
  bool operator!=(const FunctionIdentifier& other) const 
    {return *this < other || other < *this;}

  /** Return true if I am a vector */
  bool isVector() const {return algSpec().isVector();}

  /** Return true if I am a coordinate component */
  bool isCoordinateComponent() const {return algSpec().isCoordinateComponent();}

  /** Return true if I am a normal component */
  bool isNormalComponent() const {return algSpec().isNormal();}

  /** Return true if I am a scalar */
  bool isScalar() const {return algSpec().isScalar();}

private:

  /** Generate a unique ID */
  static int nextID() {static int id=0; id++; return id;}

  int dofID_;

  AlgebraSpecifier algSpec_;
};

/** \relates FunctionIdentifier */
FunctionIdentifier makeFuncID(int tensorOrder);

}


namespace std
{
/** \relates FunctionIdentifier */
ostream& operator<<(std::ostream& os, const Sundance::FunctionIdentifier& fid);
}

#endif
