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

#ifndef QUASI_RANGE_SPACE_STEP_STRATEGY_H
#define QUASI_RANGE_SPACE_STEP_STRATEGY_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Abstract interface for a strategy object that will compute a step that will
 * approximalty solve a range space subproblem {abstract}.
 */
class QuasiRangeSpaceStep_Strategy {
public:

  /** \brief . */
  virtual ~QuasiRangeSpaceStep_Strategy() {}

  /** \brief Compute a step that will approximatly solve a range-space subproblem.
   *
   * This function will compute a step <tt>v</tt> that will approximatly satisfy:
   *
   * <tt>||Gc_k'*v + c(xo)|| < ||c(xo)||</tt>
   *
   * The above norm ||.|| could be any valid norm and the implementation is free to
   * define what descent means any way it would like.  It is assumed that this
   * step will be computed by using the <tt>Gc_k</tt> but other implementations are
   * possible.  Any information being used in the algorithm can be used
   * to compute this step in a reasonable way.  Note that the 
   * inequalities do not have to (and should not in most cases) be considered
   * in this computation.  Note that whatever means is used to compute <tt>v</tt>
   * that it better give a descent direction for <tt>||c(x)||</tt> but there is no guarantee
   * for this if <tt>||xo - x_k||</tt> is large since <tt>Gc_k</tt> may not accurately approximate
   * <tt>Gc(xo)</tt>.
   *
   * @param out      [out] Output stream journal data is written to.
   * @param olevel   [in] Output level for printing to <tt>out</tt>.
   * @param algo     [in/out] The NLPAlgo object.  This object can be queryed for information.
   * @param s        [in/out] NLPAlgoState object.  May be queried or modified if needed.
   * @param xo       [in] Base point vector (size n) xo.
   * @param c_xo     [out] Constraints residual c(xo).
   * @param v        [out] Computed step vector (size n).  Must not be NULL.
   *
   * @return Returns true if a step could be found and false otherwise.
   */
   virtual bool solve_quasi_range_space_step(
    std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
    ,const Vector& xo, const Vector& c_xo, VectorMutable* v
      ) = 0;

  /** \brief This function will print a description of the computations and logic used.
   */
  virtual void print_step( std::ostream& out, const std::string& leading_str ) const = 0;

}; // end class QuasiRangeSpaceStep_Strategy

} // end namespace MoochoPack

#endif // QUASI_RANGE_SPACE_STEP_STRATEGY_H
