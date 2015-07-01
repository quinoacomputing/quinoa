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

#ifndef EVAL_NEW_POINT_STD_STEP_SET_OPTIONS_H
#define EVAL_NEW_POINT_STD_STEP_SET_OPTIONS_H

#include "MoochoPack_EvalNewPointStd_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for EvalNewPointStd_Step from an \c OptionsFromStream object.
 *
 * The default options group name is EvalNewPointStd.
 *
 * The options group is:
 *
 \verbatim
    options_group EvalNewPointStd {
        fd_deriv_testing = FD_DEFAULT;
        decomp_sys_teting = DST_DEFAULT;
        decomp_sys_teting_print_level = DSPL_USE_GLOBAL;
    }
 \verbatim
 *
 * <ul>
 * <li> <b>fd_deriv_testing</b>: Determines if finite differerece testing of the 
 *      derivatives of the Gc and Gf.  See the class \c EvalNewPointStd_Step
 *      and its printed algorithm for more details.
 *      <ul>
 *      <li> <b>FD_DEFAULT</b>: The global flag check_results determines
 *           if the tests are performed.
 *      <li> <b>FD_TEST</b>: The tests are performed reguardless the
 *           value of check_results
 *      <li> <b>FD_NO_TEST</b>: The tests are not performed reguardless the
 *           value of check_results
 *      </ul>
 * <li> <b>decomp_sys_testing</b>: Determines if the range/null decomposition of
 *      Gc and Gh is performed.  See the class \c EvalNewPointStd_Step
 *      and its printed algorithm for more details.
 *      <ul>
 *      <li> <b>DST_DEFAULT</b>: The global flag check_results determines
 *           if the tests are performed.
 *      <li> <b>DST_TEST</b>: The tests are performed reguardless the
 *           value of check_results
 *      <li> <b>DST_NO_TEST</b>: The tests are not performed reguardless the
 *           value of check_results
 *      </ul>
 * </ul>
 */
class EvalNewPointStd_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      EvalNewPointStd_Step >
{
public:

  /** \brief . */
  EvalNewPointStd_StepSetOptions(
      EvalNewPointStd_Step* target = 0
    , const char opt_grp_name[] = "EvalNewPointStd" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class EvalNewPointStd_StepSetOptions

}	// end namespace MoochoPack

#endif	// EVAL_NEW_POINT_STD_STEP_SET_OPTIONS_H
