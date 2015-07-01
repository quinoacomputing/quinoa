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

#ifndef LINE_SEARCH_WATCH_DOG_STEP_H
#define LINE_SEARCH_WATCH_DOG_STEP_H

#include "rSQPAlgo_StepBaseClasses.h"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "ConstrainedOptPack_DirectLineSearch_Strategy.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "MiStandardAggregationMacros.h"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Implements watchdog line search.
  *
  * The watchdog proceedure will only be considered when :
  * opt_kkt_err_k <= opt_kkt_err_threshold && feas_kkt_err_k <= feas_kkt_err_threshold
  * (see step listing).  The default behavior is to never use the watchdog procedure.
  */
class LineSearchWatchDog_Step : public LineSearch_Step {
public:

  /// <<std comp>> members for direct_line_search
  STANDARD_COMPOSITION_MEMBERS(DirectLineSearch_Strategy,direct_line_search);

  /// <<std comp>> members for merit_func
  STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func);

  /** \brief <<std member comp>> members for the armijo fractional reduction parameter.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, eta );

  /** \brief <<std member comp>> members for the threshold for opt_kkt_err before
    * the watchdog procedure should kick-in.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_kkt_err_threshold );
  
  /** \brief <<std member comp>> members for the threshold for feas_kkt_err before
    * the watchdog procedure should kick-in.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_kkt_err_threshold );

  /** \brief . */
  LineSearchWatchDog_Step(
        const direct_line_search_ptr_t&	direct_line_search		= 0
      , const merit_func_ptr_t&			merit_func				= 0
      , value_type						eta						= 1e-4
      , value_type						opt_kkt_err_threshold 	= 1e-1
      , value_type						feas_kkt_err_threshold 	= 1e-3
      );

  // ////////////////////
  // Overridden

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

private:
  int					watch_k_;	// >= 0 means that we are using the watchdog.
  DVector				xo_;
  value_type			fo_;
  value_type			nrm_co_;
  DVector				do_;
  value_type			phio_;
  value_type			Dphio_;
  value_type			phiop1_;
  
};	// end class LineSearchWatchDog_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_WATCH_DOG_STEP_H
