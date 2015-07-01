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

#include "MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 4;

  const char options_group_name[] = "MeritFuncPenaltyParamUpdate";

  enum local_EOptions {
      SMALL_MU,
      MIN_MU_RATIO,
      MULT_FACTOR,
      KKT_NEAR_SOL
  };

  const char* local_SOptions[local_num_options]	= {
      "small_mu",
      "min_mu_ratio",
      "mult_factor",
      "kkt_near_sol"
  };

}

namespace MoochoPack {

MeritFunc_PenaltyParamUpdate_AddedStepSetOptions::MeritFunc_PenaltyParamUpdate_AddedStepSetOptions(
      MeritFunc_PenaltyParamUpdate_AddedStep* target )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        options_group_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      MeritFunc_PenaltyParamUpdate_AddedStep >( target )
{}

void MeritFunc_PenaltyParamUpdate_AddedStepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  switch( (local_EOptions)option_num ) {
    case SMALL_MU: {
      target().small_mu( std::atof( option_value.c_str() ) );
      break;
    }
    case MIN_MU_RATIO: {
      target().min_mu_ratio( std::atof( option_value.c_str() ) );
      break;
    }
    case MULT_FACTOR: {
      target().mult_factor( std::atof( option_value.c_str() ) );
      break;
    }
    case KKT_NEAR_SOL: {
      target().kkt_near_sol( std::atof( option_value.c_str() ) );
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
