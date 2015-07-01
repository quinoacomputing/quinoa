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

#ifndef UPDATE_REDUCED_SIGMA_STEP_H
#define UPDATE_REDUCED_SIGMA_STEP_H

#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Standard class for updating the reduced sigma 
 *   for interior point optimization
 *
 *
 */

class UpdateReducedSigma_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
  {
  public:

    enum e_update_methods
      {
      ALWAYS_EXPLICIT,
      BFGS_PRIMAL,
      BFGS_DUAL_NO_CORRECTION,
      BFGS_DUAL_EXPLICIT_CORRECTION,
      BFGS_DUAL_SCALING_CORRECTION
      };

    /** \brief update method for the reduced sigma term
     *  update_method = always_explicit;
     *	update_method = BFGS_primal;
     *	update_method = BFGS_dual_no_correction;
     *	update_method = BFGS_dual_explicit_correction; *** (default)
     *	update_method = BFGS_dual_scaling_correction;
     *** These options determine exactly how the reduced sigma
     *** term will be updated. 
     ***
     *** always_explicit 				: the full Z_kT*Sigma*Zk at each step (expensive)
     *** BFGS_primal 					: a BFGS update of mu*X^-2 (exact at solution)
     *** BFGS_dual_no_correction 		: update with Z_kT*Sigma*Z_k*pz 
     ***										(no correction when mu changes)
     *** BFGS_dual_explicit_correction 	: same as above
     ***										(do an explicit calculation when mu changes)
     *** BFGS_dual_scaling_correction    : same as above
     ***										(scale by mu_kp1/mu_k when mu changes)
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( e_update_methods, update_method);

    /** @name Overridden from AlgorithmStep */
    //@{
    /** \brief . */
    bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
           , poss_type assoc_step_poss);
    
    
    void print_step( const IterationPack::Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
             , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
    //@}

    /** Constructor.
     */
    UpdateReducedSigma_Step(
//		  const e_update_methods update_method = BFGS_DUAL_EXPLICIT_CORRECTION
      const e_update_methods update_method = ALWAYS_EXPLICIT // For now only!
      );
    //@}

  private:
    void FormReducedSigmaExplicitly(NLPAlgo& algo, IpState& s, EJournalOutputLevel olevel,  std::ostream& out);


  }; // end class EvalNewPointBarrier_Step

const char UpdateReducedSigma_opt_grp_name[] = "UpdateReducedSigma";
class UpdateReducedSigma_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode,
    public OptionsFromStreamPack::SetOptionsToTargetBase< UpdateReducedSigma_Step >
  {
  public:
    UpdateReducedSigma_StepSetOptions(
      UpdateReducedSigma_Step* target = 0,
      const char opt_grp_name[] = UpdateReducedSigma_opt_grp_name );

  protected:

    /// Overridden from SetOptionsFromStreamNode
    void setOption( int option_num, const std::string& option_value );
  
  };	// end class UpdateReducedSigma_StepSetOptions


}  // end namespace MoochoPack

#endif // UPDATE_REDUCED_SIGMA_STEP_H
