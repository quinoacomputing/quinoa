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

#include "Moocho_ConfigDefs.hpp"

#ifdef HAVE_MOOCHO_MA28

#include "AbstractLinAlgPack_DirectSparseSolverMA28SetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 8;

  enum local_EOptions {
    ESTIMATED_FILLIN_RATIO
      ,U
    ,GROW
    ,TOL
    ,NSRCH
    ,LBIG
    ,PRINT_MA28_OUTPUTS
    ,OUTPUT_FILE_NAME
  };

  const char* local_SOptions[local_num_options]	= {
    "estimated_fillin_ratio"
      ,"u"
    ,"grow"
    ,"tol"
    ,"nsrch"
    ,"lbig"
    ,"print_ma28_outputs"
    ,"output_file_name"
  };

}

namespace AbstractLinAlgPack {

DirectSparseSolverMA28SetOptions::DirectSparseSolverMA28SetOptions(
      DirectSparseSolverMA28* qp_solver )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        "DirectSparseSolverMA28", local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      DirectSparseSolverMA28 >( qp_solver )
{}

void DirectSparseSolverMA28SetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;

  switch( (local_EOptions)option_num ) {
    case ESTIMATED_FILLIN_RATIO: {
      target().estimated_fillin_ratio( ::atof( option_value.c_str() ) );
      break;
    }
    case U: {
      target().u( ::atof( option_value.c_str() ) );
      break;
    }
    case GROW: {
      target().grow( StringToBool( "grow", option_value.c_str() ) );
      break;
    }
    case TOL: {
      target().tol( ::atof( option_value.c_str() ) );
      break;
    }
    case NSRCH: {
      target().nsrch( ::atoi( option_value.c_str() ) );
      break;
    }
    case LBIG: {
      target().lbig( StringToBool( "lbig", option_value.c_str() ) );
      break;
    }
    case PRINT_MA28_OUTPUTS: {
      target().print_ma28_outputs( StringToBool( "grow", option_value.c_str() ) );
      break;
    }
    case OUTPUT_FILE_NAME: {
      if( option_value == "NONE" )
        target().output_file_name( "" );
      else
        target().output_file_name( option_value );
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace AbstractLinAlgPack 

#endif // HAVE_MOOCHO_MA28
