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

#ifndef EXAMPLE_NLP_BANDED_H
#define EXAMPLE_NLP_BANDED_H


#include "NLPInterfacePack_NLPSerialPreprocessExplJac.hpp"


namespace NLPInterfacePack {


/** \brief Simple scalable serial %NLP subclass.
 *
 * This example %NLP is a scalable problem where the basis of the jacobian of the
 * equality constraints is a banded (band width = bw) symmetric positive definite
 * matrix.  Both the number of dependnet and independent variables can be varied.
 *
 * To setup this %NLP, the client specifies:<ul>
 * <li> \c nD : the number of dependent variables
 * <li> \c nI : the number of independent variables (<tt>nI <= nD</tt>)
 * <li> \c bw : the band width (defined as the number of diagonals, including the
 *              center diagonal i==j that have non-zero element)
 * <li> \c mU : the number of undecomposed dependent constraints
 * <li> \c mI : the number of general constraints
 * <li> \c diag_scal : Constant scaling factor for the diagonal elements
 * <li> \c diag_vary : The scaling factor for the diagonal elements (to produce
 *                     illconditioning) between the largets and the smallest.
 * <li> \c sym_basis : True if the basis (selected by the NLP) is symmetric or not.
 * </ul>
 *
 * This %NLP is defined as:
 \verbatim

    min    f(x) = (1/2) * sum( x(i)^2, for i = 1..n )
    s.t.
           c(j) = ( ds(j)*x(j)                                 \
                    - sum( 3/(k)*x(j-k), k=1...klu(j) )        |
                    - sum( fu/(k)*x(j+k), k=1...kuu(j) )       | for j  = 1...nD
                   ) * (x(nD+q(j)) + 1)^2                      |
                   + co(j) == 0                                /

            c(nD+jU) = c(jU) + co(nD+jU)  == 0                 } for jU = 1...mU

            hl(jI) <= x(jI) - x(nD+q(jI)) <= hu(jI)            } for jI = 1...mI

            xl(i) <= x(i) <= xu(i)                             } for i  = 1...n

    where:

        n = nD + nI

        m = nD + mU

        mI = mI

        ds(j) = diag_scal * ( (diag_vary - 1)/(nD -1) * (j - 1) + 1 )

             / 3  : if sym_basis = true
        fu = |
             \ 6  : if sym_basis = false

                                 / 2  : if floor((j-1)/nI) < nD % nI
        q(j) = floor((j-1)/nI) + |
                                 \ 1  : if floor((j-1)/nI) >= nD % nI

                                                                  
                  / bw-1 : if j - bw >= 0                         \
        klu(j) =  |                                               |
                  \ j-1  : if j - bw <= 1                         |
                                                                  | for j=1...nD
                  / bw-1 : if j + bw-1 <= nD                      |
        kuu(j) =  |                                               |
                  \ nD-j : if j - bw <= 1                         /
 \endverbatim
 * In the above formuation, the sums are not computed if the upper bounds on \c k
 * are zero.  The term <tt>co(j)</tt> is an adjustable term that can be used
 * to manipulate the solution.  Note that if <tt>co(nD+jI) != 0</tt> above, then
 * the undecomposed dependent equality constraints are inconsistent with the
 * decomposed equalities and therefore the NLP is infeasible.  An infeasible NLP
 * can also be created by manipulating \c xl(i), \c xu(i), \c hl(jI), \c hu(jI)
 * and \c co(j).
 *
 * For the above NLP, the Jacobian of the decomposed equalities has Jacobian elements:
 \verbatim
    
                      /  -3/(j-i) * (x(nD+q(j)) + 1)^2         : i - klu(i) <= j < i
                      |
                      |  ds(j) * (x(nD+q(j)) + 1)^2            : i == j
   d(c(j))/d(x(i)) =  |
                      |  -3/(i-j) * (x(nD+q(j)) + 1)^2         : i < j <= i + kuu(i)
                      |
                      |  2 * (c(j) - co(j)) / (x(nD+q(j)) + 1) : i == nD + q
                      | 
                      \  0                                     : otherwise
                      
                      , for j = 1...nD, i = 1...nD+nI
 \endverbatim
 * The above definition shows that for the independent variables, the Jacobian
 * elements are written in terms of the constraint \c c(j).  This fact
 * is exploited in the computational routines when <tt>this->multi_calc() == true</tt>.
 *
 * For <tt>nD == 7, nI == 2, bw = 2</tt> with <tt>floor(nD/nI) = 3</tt> and
 * <tt>nD \% nI = 1</tt>, the Jacobian <tt>Gc'</tt> looks like:
 \verbatim

   1 | x  x                 x    |
   2 | x  x  x              x    |
   3 |    x  x  x           x    |
   4 |       x  x  x        x    |
   5 |          x  x  x        x |
   6 |             x  x  x     x |
   7 |                x  x     x |
       -  -  -  -  -  -  -  -  -
       1  2  3  4  5  6  7  8  9
 \endverbatim
 * 
 * ToDo: Finish documentation!
 */
class ExampleNLPBanded
  : public NLPSerialPreprocessExplJac
{
public:

  /** @name Constructors / initializers */
  //@{

  /** \brief Constructor.
   *
   * ToDo: Finish documentation!
   */
  ExampleNLPBanded(
    size_type     nD
    ,size_type    nI
    ,size_type    bw                = 1
    ,size_type    mU                = 0
    ,size_type    mI                = 0
    ,value_type   xo                = 0.1
    ,value_type   xDl               = -NLP::infinite_bound()
    ,value_type   xDu               = +NLP::infinite_bound()
    ,value_type   xIl               = -NLP::infinite_bound()
    ,value_type   xIu               = +NLP::infinite_bound()
    ,value_type   hl                = -NLP::infinite_bound()
    ,value_type   hu                = +NLP::infinite_bound()
    ,bool         nlp_selects_basis = false
    ,value_type   diag_scal         = 10.0
    ,value_type   diag_vary         = 1.0
    ,bool         sym_basis         = false
    ,value_type   f_offset          = 0.0
    ,value_type   co                = 0.0
    ,bool         ignore_constraints = false
    );

  //@}

  // Todo: Add methods to manipulate bounds and co ...

  /** @name Access */
  //@{
  
  // ToDo: Add these ...

  //@}

  /** @name Overridden public members from NLP */
  //@{

  /** \brief . */
  void initialize(bool test_setup);
  /** \brief . */
  bool is_initialized() const;
  /** \brief . */
  value_type max_var_bounds_viol() const;

  //@}

  /** @name Overridden from NLPVarReductPerm */
  //@{
  
  /** \brief . */
  bool nlp_selects_basis() const;

  //@}

protected:

  /** @name Overridden protected methods from NLPSerialPreprocess */
  //@{

  /** \brief . */
  bool imp_nlp_has_changed() const;
  /** \brief . */
  size_type imp_n_orig() const;
  /** \brief . */
  size_type imp_m_orig() const;
  /** \brief . */
  size_type imp_mI_orig() const;
  /** \brief . */
  const DVectorSlice imp_xinit_orig() const;
  /** \brief . */
  bool imp_has_var_bounds() const;
  /** \brief . */
  const DVectorSlice imp_xl_orig() const;
  /** \brief . */
  const DVectorSlice imp_xu_orig() const;
  /** \brief . */
  const DVectorSlice imp_hl_orig() const;
  /** \brief . */
  const DVectorSlice imp_hu_orig() const;
  /** \brief . */
  void imp_calc_f_orig(
    const DVectorSlice            &x_full
    ,bool                        newx
    ,const ZeroOrderInfoSerial   &zero_order_info
    ) const;
  /** \brief . */
  void imp_calc_c_orig(
    const DVectorSlice            &x_full
    ,bool                        newx
    ,const ZeroOrderInfoSerial   &zero_order_info
    ) const;
  /** \brief . */
  void imp_calc_h_orig(
    const DVectorSlice            &x_full
    ,bool                        newx
    ,const ZeroOrderInfoSerial   &zero_order_info
    ) const;
  /** \brief . */
  void imp_calc_Gf_orig(
    const DVectorSlice            &x_full
    ,bool                        newx
    ,const ObjGradInfoSerial     &obj_grad_info
    ) const;
  /** \brief . */
  bool imp_get_next_basis(
    IVector      *var_perm_full
    ,IVector     *equ_perm_full
    ,size_type   *rank_full
    ,size_type   *rank
    );
  /** \brief . */
  void imp_report_orig_final_solution(
    const DVectorSlice      &x_orig
    ,const DVectorSlice     *lambda_orig
    ,const DVectorSlice     *lambdaI_orig
    ,const DVectorSlice     *nu_orig
    ,bool                  is_optimal
    );

  //@}
  
  /** @name Overridden protected methods from NLPSerialPreprocessExplJac */
  //@{

  /** \brief . */
  size_type imp_Gc_nz_orig() const;
  /** \brief . */
  size_type imp_Gh_nz_orig() const;
  /** \brief . */
  void imp_calc_Gc_orig(
    const DVectorSlice& x_full, bool newx
    , const FirstOrderExplInfo& first_order_expl_info
    ) const;
  /** \brief . */
  void imp_calc_Gh_orig(
    const DVectorSlice& x_full, bool newx
    , const FirstOrderExplInfo& first_order_expl_info
    ) const;

  //@}

private:

  // /////////////////////////////////////////
  // Private types

  // /////////////////////////////////////////
  // Private data members

  bool         is_initialized_;

  bool         nlp_selects_basis_;
  bool         basis_selection_was_given_;

  bool         has_var_bounds_;

  value_type   f_offset_;

  size_type    nD_;
  size_type    nI_;
  size_type    bw_;
  size_type    mU_;
  size_type    mI_;

  bool         ignore_constraints_;

  size_type    Gc_orig_nz_;
  size_type    Gh_orig_nz_;

  DVector      xinit_orig_;
  DVector      xl_orig_;
  DVector      xu_orig_;
  DVector      hl_orig_;
  DVector      hu_orig_;

  DVector      co_orig_;

  mutable bool c_orig_updated_;

  value_type   diag_scal_;
  value_type   diag_vary_;
  value_type   fu_;

  // /////////////////////////////////////////
  // Private member functions

  /** \brief . */
  void assert_is_initialized() const;

  /** \brief . */
  void inform_new_point(bool newx) const;

  // Not defined and not to be called
  ExampleNLPBanded();
  ExampleNLPBanded(const ExampleNLPBanded&);
  ExampleNLPBanded& operator=(const ExampleNLPBanded&);


};	// end class ExampleNLPBanded


}	// end namespace NLPInterfacePack


#endif	// EXAMPLE_NLP_BANDED_H
