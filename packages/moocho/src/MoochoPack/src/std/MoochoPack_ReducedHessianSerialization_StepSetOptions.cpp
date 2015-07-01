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

#include "MoochoPack_ReducedHessianSerialization_StepSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_Assert.hpp"

// Define the options
namespace {

const int local_num_options = 2;

enum local_EOptions {
  REDUCED_HESSIAN_INPUT_FILE_NAME
  ,REDUCED_HESSIAN_OUTPUT_FILE_NAME
};

const char* local_SOptions[local_num_options]	= {
  "reduced_hessian_input_file_name"
  ,"reduced_hessian_output_file_name"
};

inline
std::string remove_quotes( const std::string &option_name, const std::string &str )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    str[0]!='\"' || str[str.length()-1]!='\"', std::logic_error
    ,"Error, the option \'" << option_name << "\' must have a single set of quotes around it!"
    );
  if(str.length()==2)
    return std::string("");
  return str.substr(1,str.length()-2);
}

} // namespace

namespace MoochoPack {

ReducedHessianSerialization_StepSetOptions::ReducedHessianSerialization_StepSetOptions(
        ReducedHessianSerialization_Step* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      ReducedHessianSerialization_Step >( target )
{}

void ReducedHessianSerialization_StepSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;

  typedef ReducedHessianSerialization_Step target_t;
  switch( (local_EOptions)option_num ) {
    case REDUCED_HESSIAN_INPUT_FILE_NAME : {
      target().reduced_hessian_input_file_name(remove_quotes("reduced_hessian_input_file_name",option_value));
      break;
    }
    case REDUCED_HESSIAN_OUTPUT_FILE_NAME : {
      target().reduced_hessian_output_file_name(remove_quotes("reduced_hessian_output_file_name",option_value));
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace MoochoPack 
