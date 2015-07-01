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

#include <assert.h>
#include <math.h>

#include "AbstractLinAlgPack_VectorSpaceTesterSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 6;

  enum local_EOptions {
    PRINT_ALL_TESTS
    ,PRINT_VECTORS
    ,TEUCHOS_TEST_FOR_EXCEPTION
    ,NUM_RANDOM_TESTS
      ,WARNING_TOL
      ,ERROR_TOL
  };

  const char* local_SOptions[local_num_options]	= {
    "print_all_tests"
    ,"print_vectors"
    ,"throw_exception"
    ,"num_random_tests"
      ,"warning_tol"
      ,"error_tol"
  };

}

namespace AbstractLinAlgPack {

VectorSpaceTesterSetOptions::VectorSpaceTesterSetOptions(
        VectorSpaceTester* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      VectorSpaceTester >( target )
{}

void VectorSpaceTesterSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;
  typedef VectorSpaceTester target_t;
  switch( (local_EOptions)option_num ) {
    case PRINT_ALL_TESTS:
      target().print_all_tests(
        StringToBool( "print_all_tests", option_value.c_str() )
        );
      break;
    case PRINT_VECTORS:
      target().print_vectors(
        StringToBool( "print_vectors", option_value.c_str() )
        );
      break;
    case TEUCHOS_TEST_FOR_EXCEPTION:
      target().throw_exception(
        StringToBool( "throw_exception", option_value.c_str() )
        );
      break;
      case NUM_RANDOM_TESTS:
      target().num_random_tests(std::abs(std::atoi(option_value.c_str())));
      break;
      case WARNING_TOL:
      target().warning_tol(std::abs(std::atof(option_value.c_str())));
      break;
      case ERROR_TOL:
      target().error_tol(std::abs(std::atof(option_value.c_str())));
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace AbstractLinAlgPack
