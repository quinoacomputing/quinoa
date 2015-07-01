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

#include "MoochoPack_EvalNewPointStd_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_Assert.hpp"

// Define the options
namespace {

  const int local_num_options = 3;

  enum local_EOptions {
    FD_DERIV_TESTING
    ,DECOMP_SYS_TESTING
    ,DECOMP_SYS_TESTING_PRINT_LEVEL
  };

  const char* local_SOptions[local_num_options]	= {
    "fd_deriv_testing"
    ,"decomp_sys_testing"
    ,"decomp_sys_testing_print_level"
  };

}

namespace MoochoPack {

EvalNewPointStd_StepSetOptions::EvalNewPointStd_StepSetOptions(
        EvalNewPointStd_Step* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      EvalNewPointStd_Step >( target )
{}

void EvalNewPointStd_StepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;

  typedef EvalNewPointStd_Step target_t;
  switch( (local_EOptions)option_num ) {
      case FD_DERIV_TESTING:
    {
      const std::string &option = option_value.c_str();
      if( option == "FD_DEFAULT" )
        target().fd_deriv_testing( target_t::FD_DEFAULT );
      else if( option == "FD_TEST" )
        target().fd_deriv_testing( target_t::FD_TEST );
      else if( option == "FD_NO_TEST" )
        target().fd_deriv_testing( target_t::FD_NO_TEST );
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"Error, incorrect value for "
          "\"fd_deriv_testing\".  Only the options "
          "FD_DEFAULT, FD_TEST, and FD_NO_TEST "
          "are available" );
      break;
    }
      case DECOMP_SYS_TESTING:
    {
      const std::string &option = option_value.c_str();
      if( option == "DST_DEFAULT" )
        target().decomp_sys_testing( DecompositionSystemHandler_Strategy::DST_DEFAULT );
      else if( option == "DST_TEST" )
        target().decomp_sys_testing( DecompositionSystemHandler_Strategy::DST_TEST );
      else if( option == "DST_NO_TEST" )
        target().decomp_sys_testing( DecompositionSystemHandler_Strategy::DST_NO_TEST );
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"Error, incorrect value for "
          "\"decomp_sys_testing\".  Only the options "
          "DST_DEFAULT, DST_TEST, and DST_NO_TEST "
          "are available" );
      break;
    }
      case DECOMP_SYS_TESTING_PRINT_LEVEL:
    {
      const std::string &option = option_value.c_str();
      if( option == "DSPL_USE_GLOBAL" )
        target().decomp_sys_testing_print_level( DecompositionSystemHandler_Strategy::DSPL_USE_GLOBAL);
      else if( option == "DSPL_LEAVE_DEFAULT" )
        target().decomp_sys_testing_print_level( DecompositionSystemHandler_Strategy::DSPL_LEAVE_DEFAULT);
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"Error, incorrect value for "
          "\"decomp_sys_testing_print_level\".  Only the options "
          "DSPL_USE_GLOBAL and DSPL_LEAVE_DEFAULT are available" );
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
