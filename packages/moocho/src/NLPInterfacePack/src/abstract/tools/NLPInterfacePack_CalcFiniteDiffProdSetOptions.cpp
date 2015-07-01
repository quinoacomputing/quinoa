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

#include "NLPInterfacePack_CalcFiniteDiffProdSetOptions.hpp"
#include "Teuchos_Assert.hpp"

// Define the options
namespace {

  const int local_num_options = 6;

  enum local_EOptions {
    FD_METHOD_ORDER
    ,FD_STEP_SELECT
    ,FD_STEP_SIZE
    ,FD_STEP_SIZE_MIN
    ,FD_STEP_SIZE_F
    ,FD_STEP_SIZE_C
  };

  const char* local_SOptions[local_num_options]	= {
    "fd_method_order"
    ,"fd_step_select"
    ,"fd_step_size"
    ,"fd_step_size_min"
    ,"fd_step_size_f"
    ,"fd_step_size_c"
  };

}

namespace NLPInterfacePack {

CalcFiniteDiffProdSetOptions::CalcFiniteDiffProdSetOptions(
  CalcFiniteDiffProd* target
  ,const char opt_grp_name[]
  )
  :OptionsFromStreamPack::SetOptionsFromStreamNode(opt_grp_name,local_num_options,local_SOptions)
  ,OptionsFromStreamPack::SetOptionsToTargetBase<CalcFiniteDiffProd>(target)
{}

void CalcFiniteDiffProdSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  typedef CalcFiniteDiffProd target_t;
  switch( (local_EOptions)option_num ) {
      case FD_METHOD_ORDER:
    {
      const std::string &option = option_value.c_str();
      if( option == "FD_ORDER_ONE" )
        target().fd_method_order( target_t::FD_ORDER_ONE );
      else if( option == "FD_ORDER_TWO" )
        target().fd_method_order( target_t::FD_ORDER_TWO );
      else if( option == "FD_ORDER_TWO_CENTRAL" )
        target().fd_method_order( target_t::FD_ORDER_TWO_CENTRAL );
      else if( option == "FD_ORDER_TWO_AUTO" )
        target().fd_method_order( target_t::FD_ORDER_TWO_AUTO );
      else if( option == "FD_ORDER_FOUR" )
        target().fd_method_order( target_t::FD_ORDER_FOUR );
      else if( option == "FD_ORDER_FOUR_CENTRAL" )
        target().fd_method_order( target_t::FD_ORDER_FOUR_CENTRAL );
      else if( option == "FD_ORDER_FOUR_AUTO" )
        target().fd_method_order( target_t::FD_ORDER_FOUR_AUTO );
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"CalcFiniteDiffProdSetOptions::setOption(...) : Error, incorrect value for "
          "\"fd_method_order\".  Only the options FD_ORDER_ONE, FD_ORDER_TWO, "
          "FD_ORDER_TWO_CENTRAL, FD_ORDER_TWO_AUTO, FD_ORDER_FOUR, FD_ORDER_FOUR_CENTRAL "
          "and FD_ORDER_FOUR_AUTO are available" );
      break;
    }
      case FD_STEP_SELECT:
    {
      const std::string &option = option_value.c_str();
      if( option == "FD_STEP_ABSOLUTE" )
        target().fd_step_select( target_t::FD_STEP_ABSOLUTE );
      else if( option == "FD_STEP_RELATIVE" )
        target().fd_step_select( target_t::FD_STEP_RELATIVE );
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument
          ,"CalcFiniteDiffProdSetOptions::setOption(...) : Error, incorrect value for "
          "\"fd_step_select\".  Only the options are available" );
      break;
    }
      case FD_STEP_SIZE:
      target().fd_step_size(std::atof(option_value.c_str()));
      break;
      case FD_STEP_SIZE_MIN:
      target().fd_step_size_min(std::atof(option_value.c_str()));
      break;
      case FD_STEP_SIZE_F:
      target().fd_step_size_f(std::atof(option_value.c_str()));
      break;
      case FD_STEP_SIZE_C:
      target().fd_step_size_c(std::atof(option_value.c_str()));
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace NLPInterfacePack
