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

#ifndef LINE_SEARCH_FILTER_STEP_H
#define LINE_SEARCH_FILTER_STEP_H

#include <list>

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "IterationPack_CastIQMember.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "IterationPack_AlgorithmState.hpp"

namespace MoochoPack {
 
// structure for storing filter points 
class FilterEntry 
{
public:
  FilterEntry( value_type new_f, value_type new_theta, int new_iter)
    : f(new_f), theta(new_theta), iter(new_iter) {}
    
  value_type f;
  value_type theta;
  int iter;
};

// It is critical that an std::list be used because of the way iterators are
// used!
typedef std::list< FilterEntry > Filter_T;

const std::string FILTER_IQ_STRING = "LS_FilterEntries";

/** \brief Filter line-search step class.
 *
 * Todo: Finish documentataion.
 */
class LineSearchFilter_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:
    
  /** @name Public types */
  //@{

  /** \brief . */
  static value_type F_MIN_UNBOUNDED;
  
  //@}
  
  /** @name Constructors / initializers */
  //@{
  
  /** \brief Feasibility decrease fraction.
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, gamma_theta );

  /** \brief Optimality decrease fraction.
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, gamma_f );

  /** \brief Estimate of minimum value obtainable for the objective function.
   *
   * If this value is set to F_MIN_UNBOUNDED then the default behavior
   * if gamma_f is alterned otherwise the value of gamma_f used
   * is set to gamm_f_used = gamma_f *(f_k-f_min) (see the algorithm print out).
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, f_min );

  /** \brief alpha_min linearization correction fraction
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, gamma_alpha );

  /** \brief Delta parameter for switching condition.
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, delta );
  
  /** \brief Exponent for objective in switching condition.
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, s_f );

  /** \brief Exponent for theta in switching condition.
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, s_theta );

  /** \brief Factor to evaluate theta_small
   * theta_small = theta_small_fact*max(1,theta_k)
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, theta_small_fact );

  /** \brief Maximum allowable theta value
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, theta_max );

  /** \brief Constant for Armijo condition on objective
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, eta_f );

  /** \brief Backtracking fraction for step.
   *
   * ToDo: Finish documentation.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, back_track_frac );

  /** \brief Constructor.
   */
  LineSearchFilter_Step(
    Teuchos::RCP<NLPInterfacePack::NLP> nlp
    ,const std::string         obj_iq_name      = "f"
    ,const std::string         grad_obj_iq_name = "Gf"
    ,const value_type          &gamma_theta      = 1e-5
    ,const value_type          &gamma_f          = 1e-5
    ,const value_type          &f_min            = F_MIN_UNBOUNDED
    ,const value_type          &gamma_alpha      = 5e-2
    ,const value_type          &delta            = 1e-4
    ,const value_type          &s_theta          = 1.1
    ,const value_type          &s_f              = 2.3
    ,const value_type          &theta_small_fact = 1e-4 
    ,const value_type          &theta_max        = 1e10
    ,const value_type          &eta_f            = 1e-4
    ,const value_type          &back_track_frac  = 0.5
    );

  //@}
  
  /** \brief Destructor.
   */
  ~LineSearchFilter_Step();

  //@}

  /** @name Overridden from AlgorithmStep */
  //@{
  /** \brief . */
  bool do_step( Algorithm& algo, poss_type step_poss,
                IterationPack::EDoStepType type,
                poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss,
                   IterationPack::EDoStepType type,
                   poss_type assoc_step_poss, std::ostream& out,
                   const std::string& leading_str ) const;
  //@}

private:

  // /////////////////////////////
  // Private data members

  // Private Data
  CastIQMember<Filter_T> filter_;

  /// Iteration quantity access for objective value
  CastIQMember<value_type> obj_f_;

  /// ITeration quantity access for objective gradient
  CastIQMember<VectorMutable> grad_obj_f_;

  // nlp to use for calculations
  Teuchos::RCP<NLPInterfacePack::NLP> nlp_;

  // /////////////////////////////
  // Private member functions

  // Validate input parameters - fix if possible
  bool ValidatePoint(
    const IterQuantityAccess<VectorMutable>& x,
    const IterQuantityAccess<value_type>& f,
    const IterQuantityAccess<VectorMutable>* c,
    const IterQuantityAccess<VectorMutable>* h,
    const bool throw_excpt
    ) const;

  // Check that new point is not within taboo region of filter
  bool CheckFilterAcceptability(
    const value_type f, 
    const value_type theta,
    const AlgorithmState& s
    ) const;
  
  // Check the Armijo condition on f
  bool CheckArmijo(
    const value_type Gf_t_dk, 
    const value_type alpha_k, 
    const IterQuantityAccess<value_type>& f_iq
    ) const;

  // Check if f or c has sufficient reduction
  bool CheckFractionalReduction(
    const IterQuantityAccess<value_type>& f_iq,
    const value_type gamma_f_used,
    const value_type theta_kp1, 
    const value_type theta_k
    ) const;

  // Calculate the new point given d and alpha
  void UpdatePoint(
    const VectorMutable& d,
    value_type alpha, 
    IterQuantityAccess<VectorMutable>& x,
    IterQuantityAccess<value_type>& f,
    IterQuantityAccess<VectorMutable>* c,
    IterQuantityAccess<VectorMutable>* h,
    NLP& nlp
    ) const;

  // Calculate the minimum alpha before hitting restoration phase
  value_type CalculateAlphaMin(
    const value_type gamma_f_used,
    const value_type Gf_t_dk, 
    const value_type theta_k,
    const value_type theta_small
    ) const;

  // Calculate the constraint norm
  value_type CalculateTheta_k(
    const IterQuantityAccess<VectorMutable>* c,
    const IterQuantityAccess<VectorMutable>* h,
    int k
    ) const;

  // Calculate the value of gamma_k to use
  value_type CalculateGammaFUsed(
    const IterQuantityAccess<value_type> &f,
    const value_type theta_k,
    const EJournalOutputLevel olevel,
    std::ostream &out
    ) const;

  // decide if we should switch to Armijo for objective
  bool ShouldSwitchToArmijo(
    const value_type Gf_t_dk,
    const value_type alpha_k,
    const value_type theta_k,
    const value_type theta_small
    ) const;

  // Update the filter from the last iteration
  void UpdateFilter( IterationPack::AlgorithmState& s ) const;

  // Update the filter from the last iteration and Augment it with
  // the new point
  void AugmentFilter(
    const value_type gamma_f_used,
    const value_type f_with_boundary,
    const value_type theta_with_boundary,
    IterationPack::AlgorithmState& s,
    const EJournalOutputLevel olevel,
    std::ostream &out
    ) const;
    
};	// end class LineSearchFilter_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_FILTER_STEP_H
