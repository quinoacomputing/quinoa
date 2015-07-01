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

#ifndef STRING_TO_INT_MAP_H
#define STRING_TO_INT_MAP_H

#include "OptionsFromStreamPack_OptionsFromStreamExceptions.hpp"
#include "Teuchos_Assert.hpp"

namespace OptionsFromStreamPack {

/** \brief Map a string to an enumeration.
  *
  * The purpose of this class is to simplify mapping a standard string
  * to an integer which can be interpreted as an enumeration.
  *
  * Here is an example of its use.
  \verbatim

  const int n_opt = 3;
  enum MyOptEnum {
    OPT_ONE
    ,OPT_TWO
    ,OPT_THREE
  };	// must be 0, 1,..., n_opt - 1
  const char* MyOptStrings[n_opt] = {
    "OPT_ONE
    ,"OPT_TWO"
    ,"OPT_THREE"
  }; // parallels MyOptEnum
  StringToIntMap my_enum_map( "opt_map", n_opt, NyOptStrings );
  ...
  switch( (MyEnum)my_enum_map( "OPT_ONE" ) ) {
    case OPT_ONE:
      // do stuff
    case OPT_TWO:
      // do stuff
    case OPT_THREE:
      // do stuff
    default:
      // ???
  }

  \endverbatim
  * 
  * The number of strings passed to the constructor must equal the number
  * of options in the enumeration.  If there are duplicate strings
  * (capitalization concidered) then the exception #AlreadyExists# is
  * throw.  If a string that was not passed in the
  * constructor if given to #operator()( const std::string& str )# then
  * the exception #DoesNotExist# is thrown.
  *
  * In the constructor, #name# is used in error messages in the exceptions
  * thrown to help make since out of the message.
  *
  * The default constructor is not defined and not to be called.
  */
class StringToIntMap {
public:

  /** \brief . */
  class AlreadyExists : public std::logic_error
  {public: AlreadyExists(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /** \brief . */
  class DoesNotExist : public AccessException
  {public: DoesNotExist(const std::string& what_arg) : AccessException(what_arg) {}};

  /** \brief . */
  StringToIntMap( const std::string& name, int n, const char* strings[] );

  /** \brief . */
  int operator()( const std::string& str ) const;

  /** \brief . */
  const std::string& name() const;

private:
  typedef std::map< std::string, int > map_t;	// all share implementation.
  std::string name_;
  map_t map_;

  // not defined and not to be called.
  StringToIntMap();

};	// end class StringToIntMap

// ////////////////////////////////////////////
// Inline declarations

inline
const std::string& StringToIntMap::name() const
{
  return name_;
}

}	// end namespace OptionsFromStreamPack 

#endif	// STRING_TO_INT_MAP_H
