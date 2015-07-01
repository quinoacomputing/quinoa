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

#ifndef STANDARD_AGGRAGATION_MACROS_H
#define STANDARD_AGGRAGATION_MACROS_H

#include "StandardCompositionRelationshipsPack.hpp"

/** \brief \defgroup StandardAggregationMacros_grp Macros that add <<std aggr>> members for an association.
 * \ingroup Misc_grp
 *
 * For example, if you want to include a <<std aggr>> association
 * with an object of type MyClass of the name my_object you
 * would include the macro in the public section of YourClass
 * declaration as follows:
 *
 \verbatim
  class YourClass {
  public:
    STANDARD_AGGREGATION_MEMBERS( MyClass, my_object )
  };
 \endverbatim
 *
 * Note that the macro addes the private member #TYPE* NAME_#
 * to the class declaration and therefore the member NAME_ is
 * available for direct access (in a constructor for example).
 *
 * In order to have a const only association use:
 \verbatim
  class YourClass {
  public:
    STANDARD_CONST_AGGREGATION_MEMBERS( MyClass, my_object )
  };
 \endverbatim
*/
//@{

/// Insert class members for a non-const association
#define STANDARD_AGGREGATION_MEMBERS( TYPE, NAME )						\
public:																	\
  void set_ ## NAME ( TYPE* NAME )									\
  {	NAME ## _ = NAME; }												\
  TYPE* get_ ## NAME()												\
  {	return NAME ## _; }												\
  const TYPE* get_ ## NAME() const									\
  {	return NAME ## _; }												\
  TYPE& NAME()														\
  {																	\
    return StandardCompositionRelationshipsPack::role_name(			\
      NAME ## _, false, " ## NAME ## " );							\
  }																	\
  const TYPE& NAME() const											\
  {																	\
    return StandardCompositionRelationshipsPack::role_name(			\
      NAME ## _, false, " ## NAME ## " );							\
  }																	\
private:																\
  TYPE* NAME ## _;													\
public:

/// Insert class members for a constant association.
#define STANDARD_CONST_AGGREGATION_MEMBERS( TYPE, NAME )				\
public:																	\
  void set_ ## NAME ( const TYPE* NAME )								\
  {	NAME ## _ = NAME; }												\
  const TYPE* get_ ## NAME() const									\
  {	return NAME ## _; }												\
  const TYPE& NAME() const											\
  {																	\
    return StandardCompositionRelationshipsPack::const_role_name(	\
      NAME ## _, false, " ## NAME ## " );							\
  }																	\
private:																\
  const TYPE* NAME ## _;												\
public:
  
#endif	// STANDARD_AGGRAGATION_MACROS_H
