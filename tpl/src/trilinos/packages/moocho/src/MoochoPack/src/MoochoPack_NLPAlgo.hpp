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

#ifndef RSQP_ALGO_H
#define RSQP_ALGO_H

#include "MoochoPack_NLPAlgoInterface.hpp"
#include "MoochoPack_NLPAlgoContainer.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "IterationPack_Algorithm.hpp"
#include "StandardAggregationMacros.hpp"

namespace MoochoPack {

/** \brief rSQP Algorithm control class.
  */
class NLPAlgo
  : public NLPAlgoInterface
  , public IterationPack::Algorithm
{
public:

  /** @name Public Types */
  //@{

  //@}

  /// Constructs with no step, added_step, pre_step, post_step, state, or decomp_sys objects added.
  NLPAlgo();

  /// <<std aggr>> members for algo_cntr
  STANDARD_AGGREGATION_MEMBERS( NLPAlgoContainer, algo_cntr )

  /// <<std aggr>> members for nlp
  STANDARD_AGGREGATION_MEMBERS( NLP, nlp )

  /** \brief . */
  NLPAlgoState& rsqp_state()
  {	return dynamic_cast<NLPAlgoState&>(state()); }

  /** \brief . */
  const NLPAlgoState& rsqp_state() const
  {	return dynamic_cast<const NLPAlgoState&>(state()); }

  /** \brief . */
  void do_step_first(Algorithm::poss_type first_step_poss)
  {	first_step_poss_ = first_step_poss; }

  /** @name Overridden form rSQPAlgoInteface */
  //@{	
  
  /** \brief . */
  const NLPAlgoState& retrieve_state() const;

  /** \brief This is the main control function for the rSQP algorithm.
    *
    * This function basically just calls Algorithm::do_algorithm(...).
    */
  NLPSolverClientInterface::EFindMinReturn dispatch();

  /** \brief . */
  void interface_print_algorithm(std::ostream& out) const;
  /** \brief . */
  void interface_set_algo_timing( bool algo_timing );
  /** \brief . */
  bool interface_algo_timing() const;
  /** \brief . */
  void interface_print_algorithm_times( std::ostream& out ) const;

  //@}

  /// overridden from Algorihth.

  /** \brief . */
  void print_algorithm(std::ostream& out) const;

protected:

  // First step to execute
  Algorithm::poss_type first_step_poss_;

};	// end class NLPAlgo

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_H
