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

#ifndef DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_H
#define DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_H

#include "ConstrainedOptPack_DirectLineSearch_Strategy.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Performs a line search using the Armijo condition and
  * uses quadratic interpolation to select each new alpha.
  */
class DirectLineSearchArmQuad_Strategy : public DirectLineSearch_Strategy {
public:

  /// Set the Armijo cord test fractional reduction parameter.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, eta );
  
  /// The minimum fraction that alpha is reduced for each line search iteration.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, min_frac );

  /// The maximum fraction that alpha is reduced for each line search iteration.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_frac );

  /** \brief Deterimine if the line search iterations are maxed out or not.
   * 
   * This option is really only used for debugging and requires
   * changing the other parameters to make it useful.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, max_out_iter );

  /// Constructs with default settings.
  DirectLineSearchArmQuad_Strategy(
    int           max_iter       = 20
    ,value_type   eta            = 1.0e-4
    ,value_type   min_frac       = 0.1
    ,value_type   max_frac       = 0.5
    ,bool         max_out_iter   = false
    );

  /** @name Overridden from DirectLineSearch_Strategy */
  //@{

  /** \brief . */
  void set_max_iter(int max_iter);
  /** \brief . */
  int max_iter() const;
  /** \brief . */
  int num_iterations() const;
  /** \brief Performs the following line search:<br>
   *
   \verbatim

   num_iter = 0;
   while( phi.value(alpha_k) > phi_k + eta * alpha_k * phi.deriv() ) 
   {
      if(num_iter >= max_iter) return true;
      num_iter = num_iter + 1;
      alpha_k = [ min_frac * alpha_k <= quadradic interpolation for alpha	<= max_frac * alpha_k ];
   }
   return true;<br>
   \endverbatim
   * If the maximum number of iterations is exceeded then false will be returned.
   *
   * The default values for the adjustable parameters (from D&S A6.3.1)
   * are:<br>
   * max_iter = 20<br>
   * eta = 1.0e-4<br>
   * min_frac = 0.1<br>
   * max_frac = 0.5<br>
   */
  bool do_line_search(
    const MeritFuncCalc1D   &phi
    ,value_type             phi_k
    ,value_type             *alpha_k
    ,value_type             *phi_kp1
    ,std::ostream           *out
    );

  /** \brief . */
  void print_algorithm(std::ostream& out, const std::string& leading_str) const;

  //@}

private:
  int	max_iter_;
  int	num_iter_;	// stores the number of interations

  // Throw an exception if the parameters are not in a proper range.
  void validate_parameters() const;

};	// end class DirectLineSearchArmQuad_Strategy

}	// end namespace ConstrainedOptPack

#endif	// DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_H
