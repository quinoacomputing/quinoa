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

#ifndef DECOMPOSITION_SYSTEM_TESTER_SET_OPTIONS_H
#define DECOMPOSITION_SYSTEM_TESTER_SET_OPTIONS_H

#include "ConstrainedOptPack_DecompositionSystemTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace ConstrainedOptPack {

/** \brief Set options for DecompositionSystemTester from an
  * OptionsFromStream object.
  *
  * The default options group name is DecompositionSystemTester.
  *
  * The options group is:
  *
  \verbatim

    options_group DecompositionSystemTester {
        print_tests = PRINT_NONE;
    *    print_tests = PRINT_BASIC;
    *    print_tests = PRINT_MORE;
    *    print_tests = PRINT_ALL;
    *    dump_all = true;
        dump_all = false;
        throw_exception = true;
    *    throw_exception = false;
        num_random_tests = 1;
        mult_warning_tol = 1e-14;
        mult_error_tol   = 1e-10;
        solve_warning_tol = 1e-14;
        solve_error_tol   = 1e-10;
    }
  \endverbatim
  */
class DecompositionSystemTesterSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      DecompositionSystemTester >
{
public:

  /** \brief . */
  DecompositionSystemTesterSetOptions(
      DecompositionSystemTester* target = 0
    , const char opt_grp_name[] = "DecompositionSystemTester" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class DecompositionSystemTesterSetOptions

}	// end namespace ConstrainedOptPack

#endif	// DECOMPOSITION_SYSTEM_TESTER_SET_OPTIONS_H
