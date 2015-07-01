#if 0

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

#include "MoochoPack_ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 4;

  enum local_EOptions {
    MIN_NUM_UPDATES_PROJ_START
    ,MAX_NUM_UPDATES_PROJ_START
    ,NUM_SUPERBASICS_SWITCH_DENSE
    ,NUM_ADD_RECENT_UPDATES
  };

  const char* local_SOptions[local_num_options]	= {
    "min_num_updates_proj_start"
    ,"max_num_updates_proj_start"
    ,"num_superbasics_switch_dense"
    ,"num_add_recent_updates"
  };

}

namespace MoochoPack {

ReducedHessianSecantUpdateLPBFGS_StrategySetOptions::ReducedHessianSecantUpdateLPBFGS_StrategySetOptions(
  ReducedHessianSecantUpdateLPBFGS_Strategy* target
  , const char opt_grp_name[] )
  : OptionsFromStreamPack::SetOptionsFromStreamNode(
    opt_grp_name, local_num_options, local_SOptions )
  , OptionsFromStreamPack::SetOptionsToTargetBase< ReducedHessianSecantUpdateLPBFGS_Strategy >( target )
{}

void ReducedHessianSecantUpdateLPBFGS_StrategySetOptions::setOption(
  int option_num, const std::string& option_value )
{
  switch( (local_EOptions)option_num ) {
    case MIN_NUM_UPDATES_PROJ_START: {
      target().min_num_updates_proj_start( ::abs( ::atoi( option_value.c_str() ) ) );
      break;
    }
    case MAX_NUM_UPDATES_PROJ_START: {
      target().max_num_updates_proj_start( ::abs( ::atoi( option_value.c_str() ) ) );
      break;
    }
    case NUM_SUPERBASICS_SWITCH_DENSE: {
      target().num_superbasics_switch_dense( ::abs( ::atoi( option_value.c_str() ) ) );
      break;
    }
    case NUM_ADD_RECENT_UPDATES: {
      target().num_add_recent_updates( ::abs( ::atoi( option_value.c_str() ) ) );
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 

#endif // 0
