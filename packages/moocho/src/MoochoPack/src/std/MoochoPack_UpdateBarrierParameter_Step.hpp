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

#ifndef UPDATE_BARRIER_PARAMETER_STEP_H
#define UPDATE_BARRIER_PARAMETER_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Barrier Parameter (mu) Update 
 *
 * This class updates barrier_parameter & e_tol for next iteration
 *
 */

class UpdateBarrierParameter_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
  {
  public:

    /** @name Constructors / initializers */
    //@{

    /** \brief Initial barrier parameter
     *
     * mu_kp1 = min(tau_mu*mu_k,mu_k^theta_mu)
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, init_barrier_parameter );

    /** \brief barrier_parameter decrease fraction (linear decrease)
     *
     * mu_kp1 = min(tau_mu*mu_k,mu_k^theta_mu)
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, tau_mu );

    /** \brief barrier_parameter decrease power (for superlinear decrease)
     *
     * mu_kp1 = min(tau_mu*mu_k,mu_k^theta_mu)
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, theta_mu );

    /** \brief error tolerance fraction
     *
     * e_tol = min( e_tol_max, tau_epsilon*min(mu_k,mu_k^theta_epsilon))
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, tau_epsilon );

    /** \brief error tolerance power
     *
     * e_tol = min( e_tol_max, tau_epsilon*min(mu_k,mu_k^theta_epsilon))
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, theta_epsilon );

    /** \brief maximum error tolerance
     *
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, e_tol_max );

    /** \brief Constructor.
     */
    UpdateBarrierParameter_Step(
      const value_type init_barrier_parameter = 0.1,
      const value_type tau_mu = 0.2,
      const value_type theta_mu = 1.5,
      const value_type tau_epsilon = 10,
      const value_type theta_epsilon = 1.1,
      const value_type e_tol_max = 1000
      );
    //@}

    /** @name Overridden from AlgorithmStep */
    //@{
    /** \brief . */
    bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
           , poss_type assoc_step_poss);
    
    
    void print_step( const IterationPack::Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
             , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
    //@}

  private:
    value_type Calculate_e_tol(value_type mu);

  }; // end class UpdateBarrierParameter_Step


class UpdateBarrierParameter_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode,
    public OptionsFromStreamPack::SetOptionsToTargetBase< UpdateBarrierParameter_Step >
  {
  public:
    UpdateBarrierParameter_StepSetOptions(
      UpdateBarrierParameter_Step* target = 0,
      const char opt_grp_name[] = "UpdateBarrierParameter" );

  protected:

    /// Overridden from SetOptionsFromStreamNode
    void setOption( int option_num, const std::string& option_value );
  
  };	// end class UpdateBarrierParameter_StepSetOptions

} // end namespace MoochoPack

#endif // #if !defined UPDATE_BARRIER_PARAMETER_STEP_H
