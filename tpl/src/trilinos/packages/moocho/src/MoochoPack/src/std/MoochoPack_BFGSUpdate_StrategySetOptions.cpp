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

#include "MoochoPack_BFGSUpdate_StrategySetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 5;

  enum local_EOptions {
    RESCALE_INIT_IDENTITY
    ,USE_DAMPENING
    ,SECANT_TESTING
    ,SECANT_WARNING_TOL
    ,SECANT_ERROR_TOL
  };

  const char* local_SOptions[local_num_options]	= {
    "rescale_init_identity"
      ,"use_dampening"
    ,"secant_testing"
    ,"secant_warning_tol"
      ,"secant_error_tol"
  };

}

namespace MoochoPack {

BFGSUpdate_StrategySetOptions::BFGSUpdate_StrategySetOptions(
        BFGSUpdate_Strategy* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      BFGSUpdate_Strategy >( target )
{}

void BFGSUpdate_StrategySetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;
  typedef BFGSUpdate_Strategy target_t;
  switch( (local_EOptions)option_num ) {
      case RESCALE_INIT_IDENTITY:
      target().rescale_init_identity(
        StringToBool( "rescale_init_identity", option_value.c_str() ));
      break;
      case USE_DAMPENING:
      target().use_dampening(
        StringToBool( "use_dampening", option_value.c_str() ));
      break;
      case SECANT_TESTING:
    {
      const std::string &option = option_value.c_str();
      if( option == "DEFAULT" )
        target().secant_testing( target_t::SECANT_TEST_DEFAULT );
      else if( option == "TEST" )
        target().secant_testing( target_t::SECANT_TEST_ALWAYS );
      else if( option == "NO_TEST" )
        target().secant_testing( target_t::SECANT_NO_TEST );
      else
        throw std::invalid_argument( "Error, incorrect value for "
          "\"secant_testing\"." );
      break;
    }
      case SECANT_WARNING_TOL:
      target().secant_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case SECANT_ERROR_TOL:
      target().secant_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
