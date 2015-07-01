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

#include "MoochoPack_TangentialStepWithInequStd_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_Assert.hpp"

// Define the options
namespace {

  const int local_num_options = 4;

  enum local_EOptions {
    WARM_START_FRAC
    ,QP_TESTING
    ,PRIMAL_FEASIBLE_POINT_ERROR
    ,DUAL_FEASIBLE_POINT_ERROR
  };

  const char* local_SOptions[local_num_options]	= {
    "warm_start_frac"
    ,"qp_testing"
    ,"primal_feasible_point_error"
    ,"dual_feasible_point_error"
  };

}

namespace MoochoPack {

TangentialStepWithInequStd_StepSetOptions::TangentialStepWithInequStd_StepSetOptions(
  TangentialStepWithInequStd_Step* target
  ,const char opt_grp_name[]
  )
  :OptionsFromStreamPack::SetOptionsFromStreamNode(
     opt_grp_name, local_num_options, local_SOptions )
  ,OptionsFromStreamPack::SetOptionsToTargetBase<
    TangentialStepWithInequStd_Step >( target )
{}

void TangentialStepWithInequStd_StepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;

  typedef TangentialStepWithInequStd_Step target_t;
  switch( (local_EOptions)option_num ) {
    case WARM_START_FRAC:
      target().warm_start_frac(std::fabs(std::atof(option_value.c_str())));
      break;
      case QP_TESTING:
    {
      const std::string &option = option_value.c_str();
      if( option == "QP_TEST_DEFAULT" )
        target().qp_testing( target_t::QP_TEST_DEFAULT );
      else if( option == "QP_TEST" )
        target().qp_testing( target_t::QP_TEST );
      else if( option == "QP_NO_TEST" )
        target().qp_testing( target_t::QP_NO_TEST );
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"Error, incorrect value \'" << option << "\' for "
          "\"qp_testing\".  Only the options "
          "QP_TEST_DEFAULT, QP_TEST, and QP_NO_TEST "
          "are available" );
      break;
    }
        case PRIMAL_FEASIBLE_POINT_ERROR:
      target().primal_feasible_point_error(StringToBool("primal_feasible_point_error",option_value.c_str()));
      break;
      case DUAL_FEASIBLE_POINT_ERROR:
        target().dual_feasible_point_error(StringToBool("dual_feasible_point_error",option_value.c_str()));
      break;
      default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
