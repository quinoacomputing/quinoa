// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef SET_OPTIONS_TO_TARGET_BASE_H
#define SET_OPTIONS_TO_TARGET_BASE_H

#include "StandardCompositionRelationshipsPack.hpp"

namespace OptionsFromStreamPack {

/** \brief Templated node class manipulating a reference to a target object
  * who will have its options set..
  */
template< class T >
class SetOptionsToTargetBase {
public:

  /// 
  SetOptionsToTargetBase( T* target = 0 )
    : target_(target)
  {}

  /** @name <<std aggr>> stereotype members for target.
    */
  //@{

  /** \brief . */
  void set_target(T* target)
  {	target_ = target; }
  /** \brief . */
  T* get_target()
  {	return target_; }
  /** \brief . */
  const T* get_target() const
  {	return target_; }
  /** \brief . */
  T& target()
  {	return StandardCompositionRelationshipsPack::role_name(target_, false, "target"); }
  /** \brief . */
  const T& target() const
  {	return StandardCompositionRelationshipsPack::role_name(target_, false, "target"); }

  //@}

private:
  T* target_;

};	// end class SetOptionsToTargetBase

}	// end namespace OptionsFromStreamPack

#endif	// SET_OPTIONS_TO_TARGET_BASE_H
