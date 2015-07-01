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

#ifndef MERIT_FUNC_PENALTY_PARAMS_UPDATE_WITH_MULT_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAMS_UPDATE_WITH_MULT_ADDED_STEP_H

#include "MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStep.hpp"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Updates a set of penalty parameters for a merit function as:
  * mu(j) = max( mu(j), |lambda_k(j)| ).
  *
  * This class assumes the merit function supports the interfaces
  * MeritFuncPenaltyParams and MeritFuncNLPDirecDeriv.
  */
class MeritFunc_PenaltyParamsUpdateWithMult_AddedStep
  : public MeritFunc_PenaltyParamUpdate_AddedStep
{
public:

  /// <<std comp>> members for merit_func
  STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func);

  /** \brief . */
  MeritFunc_PenaltyParamsUpdateWithMult_AddedStep(const merit_func_ptr_t& merit_func
    , value_type small_mu = 1e-6, value_type min_mu_ratio = 1e-8
    , value_type mult_factor = 1e-4 , value_type kkt_near_sol = 1e-1 );

  // ///////////////////////////////
  // Overridden from AlgorithmStep

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss
    , IterationPack::EDoStepType type, poss_type assoc_step_poss
    , std::ostream& out, const std::string& leading_str ) const;

  // //////////////////////////////////////////////////////////////////////
  // Overridden from MeritFunc_PenaltyParamUpdate_AddedStep

  /** \brief . */
  void small_mu( value_type small_mu );
  /** \brief . */
  value_type small_mu() const;

  /** \brief . */
  void min_mu_ratio( value_type min_mu_ratio );
  /** \brief . */
  value_type min_mu_ratio() const;

  /** \brief . */
  void mult_factor( value_type mult_factor );
  /** \brief . */
  value_type mult_factor() const;

  /** \brief . */
  void kkt_near_sol( value_type kkt_near_sol );
  /** \brief . */
  value_type kkt_near_sol() const;

private:
  bool near_solution_;
  value_type small_mu_;
  value_type min_mu_ratio_;
  value_type mult_factor_;
  value_type kkt_near_sol_;
  value_type norm_inf_mu_last_;
  
};	// end class MeritFunc_PenaltyParamsUpdateWithMult_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_PENALTY_PARAMS_UPDATE_WITH_MULT_ADDED_STEP_H
