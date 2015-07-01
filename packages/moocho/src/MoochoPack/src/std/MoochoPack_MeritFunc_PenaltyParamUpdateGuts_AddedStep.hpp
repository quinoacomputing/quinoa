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

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_GUTS_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_GUTS_ADDED_STEP_H

#include "MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStep.hpp"

namespace MoochoPack {

/** \brief Updates the penalty parameter for a merit function as:
 * mu_k = max( mu_km1, min_mu ).
 *
 * This class assumes the merit function iteration quantity supports
 * the interfaces <tt>MeritFuncPenaltyParam</tt> and <tt>MeritFuncNLPDirecDeriv</tt>.
 * 
 * min_mu is computed by subclasses.
 */
class MeritFunc_PenaltyParamUpdateGuts_AddedStep
  : public MeritFunc_PenaltyParamUpdate_AddedStep
{
public:

  /** \brief . */
  MeritFunc_PenaltyParamUpdateGuts_AddedStep(
    value_type   small_mu
    ,value_type  mult_factor
    ,value_type  kkt_near_sol
    );

  /** @name Overridden from AlgorithmStep */
  //@{

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss
    , IterationPack::EDoStepType type, poss_type assoc_step_poss
    , std::ostream& out, const std::string& leading_str ) const;

  //@}

  /** @name Overridden from MeritFunc_PenaltyParamUpdate_AddedStep */
  //@{

  /** \brief . */
  void small_mu( value_type small_mu );
  /** \brief . */
  value_type small_mu() const;
  /** \brief . */
  void mult_factor( value_type mult_factor );
  /** \brief . */
  value_type mult_factor() const;
  /** \brief . */
  void kkt_near_sol( value_type kkt_near_sol );
  /** \brief . */
  value_type kkt_near_sol() const;
  
  //@}

protected:

  /** @name Pure virtual functions to be overridden by subclasses */
  //@{

  /** \brief Override to determine the mininum value of mu the penalty parameter
   * can take on.
   *
   * @param  s       [in] The rSQP state object to get at useful information.
   * @param  min_mu  [out] If min_mu(...) returns true, then this is the
   * 					mininum value mu can take on and still have
   * 					descent in the merit function.
   * 					
   * @return	Returns true if the penalty parameter should be updated
   * 	or false if the previous mu_km1 should be used.
   */
  virtual bool min_mu( NLPAlgoState& s, value_type* min_mu ) const = 0;

  /** \brief Override to print how min_mu calculated.
   */
  virtual void print_min_mu_step( std::ostream& out
    , const std::string& leading_str ) const = 0;

  //@}

private:
  bool near_solution_;
  value_type small_mu_;
  value_type mult_factor_;
  value_type kkt_near_sol_;
  
};	// end class MeritFunc_PenaltyParamUpdateGuts_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_GUTS_ADDED_STEP_H
