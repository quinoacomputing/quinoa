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

#include "ConstrainedOptPack_DirectLineSearchArmQuad_StrategySetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 5;

  enum local_EOptions {
    SLOPE_FRAC,
    MIN_FRAC_STEP,
    MAX_FRAC_STEP,
    MAX_LS_ITER,
    MAX_OUT_LS_ITER
  };

  const char* local_SOptions[local_num_options]	= {
    "slope_frac",
    "min_frac_step",
    "max_frac_step",
    "max_ls_iter",
    "max_out_ls_iter"
  };

}

namespace ConstrainedOptPack {

DirectLineSearchArmQuad_StrategySetOptions::DirectLineSearchArmQuad_StrategySetOptions(
        DirectLineSearchArmQuad_Strategy* qp_solver
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      DirectLineSearchArmQuad_Strategy >( qp_solver )
{}

void DirectLineSearchArmQuad_StrategySetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;
  switch( (local_EOptions)option_num ) {
    case SLOPE_FRAC: {
      target().eta( std::atof( option_value.c_str() ) );
      break;
    }
    case MIN_FRAC_STEP: {
      target().min_frac( std::atof( option_value.c_str() ) );
      break;
    }
    case MAX_FRAC_STEP: {
      target().max_frac( std::atof( option_value.c_str() ) );
      break;
    }
    case MAX_LS_ITER: {
      target().set_max_iter( std::atof( option_value.c_str() ) );
      break;
    }
    case MAX_OUT_LS_ITER: {
      target().max_out_iter( StringToBool( "max_out_ls_iter", option_value.c_str() ) );
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace ConstrainedOptPack 
