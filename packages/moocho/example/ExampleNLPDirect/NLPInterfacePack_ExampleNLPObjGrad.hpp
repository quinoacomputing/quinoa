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

#ifndef EXAMPLE_NLP_OBJ_GRADIENT_H
#define EXAMPLE_NLP_OBJ_GRADIENT_H

#include "NLPInterfacePack_NLPObjGrad.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorSpaceBlocked.hpp"
#include "Teuchos_Assert.hpp"

namespace NLPInterfacePack {

/** \brief Simple example %NLP subclass to illustrate how to implement the
 * \c NLPObjGrad interface for a specialized \c NLP.
 *
 * The example %NLP we will use is a scalable problem where
 * the basis of the jacobian of the constraints is a diagonal
 * matrix (however it is not computed here).
 \verbatim

    min    f(x) = (1/2) * sum( x(i)^2, for i = 1..n )
    s.t.   c(x)(j) = x(j) * (x(m+j) -1) - 10 * x(m+j) = 0, for j = 1..m
          0.01 < x(i) < 20, for i = p...p+m

    where:
        m = n/2
        p = 1 if dep_bounded == true or m+1 if dep_bounded = false
 \endverbatim
 * This is not really a fully functional NLP in the sense that there is
 * no derivative information for the constraints.
 */
class ExampleNLPObjGrad : virtual public NLPObjGrad {
public:

  /** \brief Constructor.
   *
   * @param  vec_space  [in] Smart pointer to a vector space object that will
   *                    be used to define the spaces of dependent and independent
   *                    variables.
   * @param  xo         [in] The initial starting guess for \a x.
   * @param  has_bounds [in] If \c true, then the NLP will have bounds.  If \c false
   *                    then it will not have bounds.
   * @param  dep_bouned [in] If \c true, then the bounds will be on the dependent
   *                    variables.  If \c false, then the bounds will be on the
   *                    independent variable.  This argument is ignored if
   *                    <tt>has_bounds == false</tt>.
   */
  ExampleNLPObjGrad(
    const VectorSpace::space_ptr_t&  vec_space
    ,value_type                      xo
    ,bool                            has_bounds
    ,bool                            dep_bounded
    );

  /** @name Helper methods to be used by subclasses. */
  //@{

  /** \brief . */
  virtual Range1D var_dep() const;
  /** \brief . */
  virtual Range1D var_indep() const;

  //@}

  /** @name Overridden public members from NLP */
  //@{

  /** \brief . */
  void initialize(bool test_setup);
  /** \brief . */
  bool is_initialized() const;
  /** \brief . */
  size_type n() const;
  /** \brief . */
  size_type m() const;
  /** \brief . */
  vec_space_ptr_t space_x() const;
  /** \brief . */
  vec_space_ptr_t space_c() const;
  /** \brief . */
    size_type num_bounded_x() const;
  /** \brief . */
  void force_xinit_in_bounds(bool force_xinit_in_bounds);
  /** \brief . */
  bool force_xinit_in_bounds() const;
  /** \brief . */
  const Vector& xinit() const;
  /** \brief . */
  const Vector& xl() const;
  /** \brief . */
  const Vector& xu() const;
  /** \brief . */
  value_type max_var_bounds_viol() const;
  /** \brief . */
  void scale_f( value_type scale_f );
  /** \brief . */
  value_type scale_f() const;
  /** \brief . */
  void report_final_solution(
    const Vector&    x
    ,const Vector*   lambda
    ,const Vector*   nu
    ,bool            optimal
    );

  //@}

protected:

  /** @name Overridden protected members from NLP */
  //@{

  /** \brief . */
  void imp_calc_f(
    const Vector& x, bool newx
    ,const ZeroOrderInfo& zero_order_info) const;
  /** \brief . */
  void imp_calc_c(
    const Vector& x, bool newx
    ,const ZeroOrderInfo& zero_order_info) const;
  /// This implementation does nothing (should never be called though).
  void imp_calc_h(const Vector& x, bool newx, const ZeroOrderInfo& zero_order_info) const;

  //@}

  /** @name Overridden protected members from NLPObjGrad */
  //@{

  /** \brief . */
  void imp_calc_Gf(
    const Vector& x, bool newx
    ,const ObjGradInfo& obj_grad_info) const;

  //@}

private:

  // /////////////////////////////////////////
  // Private data members

  VectorSpace::space_ptr_t    vec_space_;       // The vector space for dependent and indepenent variables and c(x).
  VectorSpace::space_ptr_t    vec_space_comp_;  // Composite vector space for x = [ xD; xI ]
  Range1D                     var_dep_;         // Range for dependnet variables.
  Range1D                     var_indep_;       // Range for independent variables.

  bool         initialized_;            // flag for if initialized has been called.
  value_type   obj_scale_;              // default = 1.0;
  bool         has_bounds_;             // default = true
  bool         force_xinit_in_bounds_;  // default = true.

  size_type    n_;                      // Number of variables in the problem.
  VectorSpace::vec_mut_ptr_t  xinit_;   // Initial guess.
  VectorSpace::vec_mut_ptr_t  xl_;      // lower bounds.
  VectorSpace::vec_mut_ptr_t  xu_;      // upper bounds.

  // /////////////////////////////////////////
  // Private member functions

  /** \brief . */
  void assert_is_initialized() const;

};	// end class ExampleNLPObjGrad

// ///////////////////////////////////////////////
// Inline member functions

inline
void ExampleNLPObjGrad::assert_is_initialized() const
{
  typedef NLPInterfacePack::NLP NLP;
  TEUCHOS_TEST_FOR_EXCEPTION(
    !is_initialized(), NLP::UnInitialized
    ,"ExampleNLPObjGrad::assert_is_initialized() : Error, "
    "ExampleNLPObjGrad::initialize() has not been called yet." );
}

}	// end namespace NLPInterfacePack

#endif	// EXAMPLE_NLP_OBJ_GRADIENT_H
