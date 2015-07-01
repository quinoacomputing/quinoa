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

#ifndef MERIT_FUNC_CALC_1D_QUADRATIC_H
#define MERIT_FUNC_CALC_1D_QUADRATIC_H

#include "ConstrainedOptPack_MeritFuncCalc1D.hpp"
#include "ConstrainedOptPack_MeritFuncCalc.hpp"

namespace ConstrainedOptPack {

/** \brief Adds the ability to compute phi(alpha) at alpha of a given set of vectors.
  *
  * Computes <tt>phi( x = sum( alpha^k * d[k], k = 0...p-1 ) )</tt> where 
  * <tt>1 <= p <= 2</tt>.
  */
class MeritFuncCalc1DQuadratic : public MeritFuncCalc1D {
public:

  /** \brief . */
  typedef const Vector*   const_VectorWithOp_ptr;

  /** @name Constructors */
  //@{

  /** \brief The only constructor.
    *
    * Note that \c *x and \c *d gets updated as <tt>operator()(alpha)</tt> is called.
    *
    * The client must ensure that the memory pointed to by the vectors in d must not
    * be desturbed while this object is in use.  To do so may have bad side effects.
    *
    * @param  phi  [in] The merit function to use.
    * @param  p    [in] The number of vectors in \c d[].
    * @param  d    [in] Array (length \c p) of pointers to the rhs \c d[] vectors.
    * @param  x    [out] The vector that gets updated.
    *
    * Preconditions:<ul>
    * <li> <tt>1 <= p <= 3<tt>
    * <li> <tt>d[k]->space().is_compatible(x->space()), for k = 0...p-1</tt>
    * </ul>
    */
  MeritFuncCalc1DQuadratic(
    const MeritFuncCalc&      phi
    ,size_type                p
    ,const_VectorWithOp_ptr   d[]
    ,VectorMutable*     x
    );

  //@}

  /** @name Overridden from MeritFuncCalc1D */
  //@{

  /// Returns <tt>phi( x = sum( alpha^k * d[k], k = 0...p-1 ) )</tt>.
  value_type operator()(value_type alpha) const;

  /// Returns phi.deriv()
  value_type deriv() const;

  /// Calls <tt>phi->print_merit_func()</tt>.
  void print_merit_func(
    std::ostream& out, const std::string& leading_str ) const;

  //@}

private:
  const MeritFuncCalc&        phi_;
  size_type                   p_;
  const_VectorWithOp_ptr      d_[3];
  VectorMutable         *x_;

  // not defined and not to be called
  MeritFuncCalc1DQuadratic();
  MeritFuncCalc1DQuadratic& operator=( const MeritFuncCalc1DQuadratic& );

};	// end class MeritFuncCalc1DQuadratic

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_1D_QUADRATIC_H
