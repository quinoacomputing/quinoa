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

#include "MoochoPack_LineSearchFilter_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_Assert.hpp"

// Define the options
namespace {

const int local_num_options = 11;

enum local_EOptions 
  {
    GAMMA_THETA
    ,GAMMA_F
    ,F_MIN
    ,GAMMA_ALPHA
    ,DELTA
    ,S_THETA
    ,S_F
    ,THETA_SMALL_FACT
    ,THETA_MAX
    ,ETA_F
    ,BACK_TRACK_FRAC
  };

const char* local_SOptions[local_num_options] = 
  {
    "gamma_theta"
    ,"gamma_f"
    ,"f_min"
    ,"gamma_alpha"
    ,"delta"
    ,"s_theta"
    ,"s_f"
    ,"theta_small_fact"
    ,"theta_max"
    ,"eta_f"
    ,"back_track_frac"
  };

}

namespace MoochoPack {

LineSearchFilter_StepSetOptions::LineSearchFilter_StepSetOptions(
  LineSearchFilter_Step* target
  , const char opt_grp_name[] )
  :
  OptionsFromStreamPack::SetOptionsFromStreamNode(
    opt_grp_name, local_num_options, local_SOptions ),
  OptionsFromStreamPack::SetOptionsToTargetBase< LineSearchFilter_Step >( target )
  {
  }

void LineSearchFilter_StepSetOptions::setOption( 
  int option_num, const std::string& option_value )
  {
  using OptionsFromStreamPack::StringToBool;
  
  typedef LineSearchFilter_Step target_t;
  switch( (local_EOptions)option_num ) {
    case GAMMA_THETA:
      target().gamma_theta(std::atof(option_value.c_str()));
      break;
    case GAMMA_F:
      target().gamma_f(std::atof(option_value.c_str()));
      break;
    case F_MIN: {
      if( option_value == "UNBOUNDED" )
        target().f_min(target_t::F_MIN_UNBOUNDED);
      else
        target().f_min(std::atof(option_value.c_str()));
      break;
    }
    case GAMMA_ALPHA:
      target().gamma_alpha(std::atof(option_value.c_str()));
      break;
    case DELTA:
      target().delta(std::atof(option_value.c_str()));
      break;
    case S_THETA:
      target().s_theta(std::atof(option_value.c_str()));
      break;
    case S_F:
      target().s_f(std::atof(option_value.c_str()));
      break;
    case THETA_SMALL_FACT:
      target().theta_small_fact(std::atof(option_value.c_str()));
      break;
    case THETA_MAX:
      target().theta_max(std::atof(option_value.c_str()));
      break;
    case ETA_F:
      target().eta_f(std::atof(option_value.c_str()));
      break;
    case BACK_TRACK_FRAC:
      target().back_track_frac(std::atof(option_value.c_str()));
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
    }
  }

}	// end namespace MoochoPack 
