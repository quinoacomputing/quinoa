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

#ifndef CHECK_DESCENT_RANGE_SPACE_STEP_STEP_H
#define CHECK_DESCENT_RANGE_SPACE_STEP_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "NLPInterfacePack_CalcFiniteDiffProd.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Checks for descent in the decomposed equality constraints with respect to the
 * range space step <tt>Ypy</tt> using finite differences.
 *
 * This step class checks for descent in the feasibility measure <tt>q(x) = 1/2 * cd(x)'*cd(x) <: REAL</tt>
 * of the decomposed equality constraints <tt>cd(x) = c(equ_decomp)(x)</tt> with respect to the range space
 * step <tt>Ypy_k</tt>.  The gradient of this feasibility measure is:
 \verbatim

 grad(q(x),x) = grad(cd(x),x) * cd(x)
 \endverbatim
 * Therefore, we can determine if we have descent by checking <tt>grad(q(x),x)'*Ypy_k = cd(x)'*grad(cd(x),x)'*Ypy_k< 0</tt>.
 * The product <tt>grad(c(x),x)'*Ypy_k</tt> is approximated with finite differences using the class
 * <tt>MoochoPack::CalcFiniteDiffProd</tt>.
 */
class CheckDescentQuasiNormalStep_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /// Set the object that will compute the finite difference products.
  STANDARD_COMPOSITION_MEMBERS( CalcFiniteDiffProd, calc_fd_prod );

  /** \brief Constructor
   */
  CheckDescentQuasiNormalStep_Step(
    const calc_fd_prod_ptr_t&   calc_fd_prod
    );

  /** @name Overridden from AlgorithmStep */
  //@{
  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
  //@}

};	// end class CheckDescentQuasiNormalStep_Step

}	// end namespace MoochoPack 

#endif	// CHECK_DESCENT_RANGE_SPACE_STEP_STEP_H
