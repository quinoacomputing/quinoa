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

#ifndef FEASIBILITY_STEP_STRATEGY_H
#define FEASIBILITY_STEP_STRATEGY_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Abstract interface for a strategy object that will compute a step that will
 * improve feasibility (at least descent) {abstract}.
 */
class FeasibilityStep_Strategy {
public:

  /** \brief . */
  virtual ~FeasibilityStep_Strategy() {}

  /** \brief Compute a step that improves feasibility (at least locally).
   *
   * This function will compute a step <tt>w</tt> that satisfies:
   *
   * <tt>d_bounds_k.l <= xo + w - x_k <= d_bounds_k.u</tt>
   *
   * and gives descent for <tt>||c(xo + beta*w)||</tt> for at least small <tt>beta > 0</tt>.
   * This norm <tt>||.||</tt> could be any valid norm and the implementation is free to
   * define what descent means any way it would like.  Any information being used
   * in the algorithm can be used to compute this step.
   *
   * @param out     [out] Output stream journal data is written to.
   * @param olevel  [in] Output level for printing to #out#
   * @param algo    [in/out] The NLPAlgo object.  This object can be queryed for
   *                information.
   * @param s       [in/out] NLPAlgoState object.  May be queried or modified if needed.
   * @param xo      [in] Base point vector (size n) xo.
   * @param c_xo    [in] c(xo).
   * @param w       [out] Computed step vector (size n) w.  Must not be NULL.
   *
   * @return Returns true if a step that reduces feasibility subject to the bounds could
   * be found and false otherwise.
   */
   virtual bool compute_feasibility_step(
    std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
    ,const Vector& xo, const Vector& c_xo, VectorMutable* w
      ) = 0;

  /** \brief This function will print a description of the computations and logic used.
   */
  virtual void print_step( std::ostream& out, const std::string& leading_str ) const = 0;

}; // end class FeasibilityStep_Strategy

} // end namespace MoochoPack

#endif // FEASIBILITY_STEP_STRATEGY_H
