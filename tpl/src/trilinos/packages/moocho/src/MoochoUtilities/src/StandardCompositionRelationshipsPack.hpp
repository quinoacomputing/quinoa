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

#ifndef STANDARD_COMPOSITION_RELATIONSHIPS_PACK_H
#define STANDARD_COMPOSITION_RELATIONSHIPS_PACK_H

#include "Moocho_ConfigDefs.hpp"
#include "Teuchos_Assert.hpp"

namespace StandardCompositionRelationshipsPack {

/** @name <<std comp>> Stereotype Implementation Helper Package.
  *
  * This is the set of helper functions specified in the diagram
  * "Class Diagram : <<std comp>> Stereotype Implementation Helper Package".
  */
//@{

/// Thrown when the reference has not been set.
class NoRefSet : public std::logic_error
{public: NoRefSet(const std::string& what_arg) : std::logic_error(what_arg) {}};

// Throw a NoRefSet exception
inline void ThrowNoRefSet(const char func_name[], const char name[])
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,NoRefSet
    ,func_name << ": The reference for \'" << name << "\' has not been set yet"
    );
}

/// Assert that the reference is set.
template<class ContainedClass>
inline void assert_role_name_set(const ContainedClass* role_name_, const char func_name[]
  , const char name[])
{
  if(!role_name_) ThrowNoRefSet(func_name, name);
}

/** \brief . */
template<class ContainedClass>
inline void set_role_name(ContainedClass*& role_name_, bool& owns_role_name_, const char name[]
  , ContainedClass* role_name, bool owns_role_name)
{
  if(owns_role_name_ && role_name_ != role_name) delete role_name_;
  role_name_ = role_name; owns_role_name_ = owns_role_name;
}

/** \brief . */
template<class ContainedClass>
inline ContainedClass* get_role_name(ContainedClass* role_name_, bool owns_role_name_
  , const char name[])
{
  return role_name_;
}

/** \brief . */
template<class ContainedClass>
inline void set_owns_role_name(ContainedClass*& role_name_, bool& owns_role_name_
  , const char name[], bool owns_role_name)
{
  assert_role_name_set(role_name_, "set_owns_role_name()", name);
  owns_role_name_ = owns_role_name;
}

/** \brief . */
template<class ContainedClass>
inline bool owns_role_name(ContainedClass* role_name_, bool owns_role_name_, const char name[])
{
  assert_role_name_set(role_name_, "owns_role_name()", name);
  return owns_role_name_;
}

/** \brief . */
template<class ContainedClass>
inline ContainedClass& role_name(ContainedClass* role_name_, bool owns_role_name_, const char name[])
{
  assert_role_name_set(role_name_, "role_name()", name);
  return *role_name_;
}

/** \brief . */
template<class ContainedClass>
inline const ContainedClass& role_name(const ContainedClass* role_name_, bool owns_role_name_, const char name[])
{
  assert_role_name_set(role_name_, "role_name()", name);
  return *role_name_;
}

/** \brief . */
template<class ContainedClass>
inline const ContainedClass& const_role_name(const ContainedClass* role_name_, bool owns_role_name_, const char name[])
{
  assert_role_name_set(role_name_, "role_name()", name);
  return *role_name_;
}


/** \brief . */
template<class ContainedClass>
inline void destory_container_obj(ContainedClass* role_name_, bool owns_role_name_)
{
  if(owns_role_name_) delete role_name_;
}

//@}

}	// end namespace StandardCompositionRelationshipsPack

#endif // STANDARD_COMPOSITION_RELATIONSHIPS_PACK_H
