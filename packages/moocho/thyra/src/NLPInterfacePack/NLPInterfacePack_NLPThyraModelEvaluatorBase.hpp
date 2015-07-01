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

#ifndef NLPIP_NLP_THYRA_MODEL_EVALUATOR_BASE_HPP
#define NLPIP_NLP_THYRA_MODEL_EVALUATOR_BASE_HPP

#include <vector>

#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack { class VectorSpaceThyra; }

namespace NLPInterfacePack {

/** \brief Implements the base %NLP interface using a
 * <tt>Thyra::ModelEvaluator</tt> object.
 *
 * The nonlinear program is mapped as follows:

 \verbatim

    min     f(x)
    s.t.    c(x) = 0
            xL <= x <= xu

    where:

        x = [ model.x        ]
            [ model.p(p_idx) ]

        f(x) = model.g(g_idx)

        c(x) = model.f()

 \endverbatim

 * where <tt>p_idx > 0</tt> and <tt>g_idx > 0</tt> are some indexes that
 * specificy the indepenent variables and the objective function.
 *
 * In addition, the client can also override the bounds on <tt>model.x</tt> and
 * <tt>model.p(p_idx)</tt> defined in the object <tt>model</tt>.
 *
 * The current implementation of this class does not allow the use of any of
 * the auxiliary functions <tt>model.g()</tt> as undecomposed equality
 * constraints or extra general inequality constraints.  This type of
 * functionality can be added when it is needed (just ask for it).
 *
 * ToDo: Finish documentation!
 */
class NLPThyraModelEvaluatorBase : virtual public NLPObjGrad {
public:

  /** \brief Set if a trace of the model evaluations is shown or not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, showModelEvaluatorTrace );

  /** @name Overridden public members from NLP */
  //@{

  /** \brief . */
  void initialize(bool test_setup);
  /** \brief . */
  bool is_initialized() const;
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
  void set_f(value_type* f);
  /** \brief . */
  void set_c(VectorMutable* c);
  /** \brief . */
  void unset_quantities();
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

  /** @name Overridden public members from NLPObjGrad */
  //@{

  /** \brief . */
  void set_Gf(VectorMutable* Gf);

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

  //@}

  /** @name Overridden protected members from NLPObjGrad */
  //@{

  /** \brief . */
  void imp_calc_Gf(
    const Vector& x, bool newx
    ,const ObjGradInfo& obj_grad_info) const;

  //@}

protected:

  /** @name Protected functions to be used by subclasses */
  //@{

  /** Initialize to uninitialized */
  NLPThyraModelEvaluatorBase();

  /** \brief Initialize given a <tt>Thyra::ModelEvaluator</tt> and
   * a description of how to interpret it.
   *
   * @param  model    [in] NonlinearProblem that defines all of the functions and variables.
   * @param  p_idx [in] Index of the subset of parameter vectors to use as the independent
   *               variables.  If <tt>p_idx < 0</tt>, then no extra parameters are added.
   * @param  g_idx [in] Index of the subset of auxiliary response functions to use as
   *               the objective function.  Note, only the first element <tt>model.g(g_idx)(1)</tt>
   *               will be used as the objective function value.
   * @param  model_xL [in] Pointer to upper bounds for the state variables <tt>model.x</tt>.  If NULL
   *               then the default supplied in <tt>model->get_x_lower_bounds()</tt> will be used.
   * @param  model_xU [in] Pointer to upper bounds for the state variables <tt>x</tt>.  If NULL
   *               then the default supplied in <tt>model->get_x_upper_bounds()</tt> will be used.
   * @param  model_x0 [in] Pointer to initial guess for the state variables <tt>x</tt>.  If NULL
   *               the the default supplied in <tt>model->get_x_init()</tt> will be used.
   *
   * ToDo: Finish documentation!
   *
   * Todo: Add arguments for auxiliary inequalites and equalities
   */
  void initializeBase(
    const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
    );

  /** \brief Update the initial guess and bounds . */
  void updateInitialGuessAndBounds() const;

  /** \brief . */
  void assert_is_initialized() const;
  /** \brief . */
  void copy_from_model_x( const Thyra::VectorBase<value_type>* model_x, VectorMutable* x_D ) const;
  /** \brief . */
  void copy_from_model_p( const Thyra::VectorBase<value_type> *model_p, VectorMutable* x_I ) const;
  /** \brief . */
  void set_x(
    const Vector                                      &x
    ,Thyra::ModelEvaluatorBase::InArgs<value_type>    *model_inArgs_inout
    ) const;
  /** \brief . */
  void preprocessBaseInOutArgs(
    const Vector                                      &x
    ,bool                                             newx
    ,const ZeroOrderInfo                              *zero_order_info
    ,const ObjGradInfo                                *obj_grad_info
    ,const NLPFirstOrder::FirstOrderInfo              *first_order_info
    ,Thyra::ModelEvaluatorBase::InArgs<value_type>    *model_inArgs_inout
    ,Thyra::ModelEvaluatorBase::OutArgs<value_type>   *model_outArgs_inout
    ,MatrixOp*                                        *Gc_out
    ,VectorMutable*                                   *Gf_out
    ,value_type*                                      *f_out
    ,VectorMutable*                                   *c_out
    ) const;
  /** \brief . */
  void postprocessBaseOutArgs(
    Thyra::ModelEvaluatorBase::OutArgs<value_type>        *model_outArgs_inout
    ,VectorMutable                                        *Gf
    ,value_type                                           *f
    ,VectorMutable                                        *c
    ) const;
    
  //@}

//private: // ToDo: Make these private and refactor the other classes ...

  // /////////////////////////////////////////
  // Private types

  typedef Teuchos::RCP<const AbstractLinAlgPack::VectorSpaceThyra> VectorSpaceThyra_ptr_t;

  // /////////////////////////////////////////
  // Private data members

  bool                                initialized_;  // flag for if initialized has been called.
  value_type                          obj_scale_;    // default = 1.0;
  bool                                has_bounds_;   // True if has bounds
  bool                                force_xinit_in_bounds_; // default = true.
  index_type                          num_bounded_x_;
  Teuchos::RCP<Thyra::ModelEvaluator<value_type> >
                                      model_;
  int                                 p_idx_;
  int                                 g_idx_;
  bool                                DfDp_supports_op_;
  bool                                DfDp_supports_mv_;
  VectorSpace::space_ptr_t            space_x_;      // Space for the variables
  VectorSpaceThyra_ptr_t              space_c_;      // Space for the constraints
  NLPFirstOrder::mat_fcty_ptr_t       factory_Gc_;   // Factory for Gc
  NLPFirstOrder::basis_sys_ptr_t      basis_sys_;    // The basis system
  mutable bool                        x_guess_bounds_updated_;
  VectorSpace::vec_mut_ptr_t          xinit_;        // Initial guess.
  VectorSpace::vec_mut_ptr_t          xl_;           // lower bounds.
  VectorSpace::vec_mut_ptr_t          xu_;           // upper bounds.

  Teuchos::RCP<Thyra::VectorBase<value_type> >                     model_g_;

  mutable bool model_g_updated_;
  mutable bool model_Dg_updated_;

  mutable bool f_updated_;
  mutable bool c_updated_;
  mutable bool Gf_updated_;
  mutable bool Gc_updated_;

  // /////////////////////////////////////////
  // Private member functions

  /** \brief . */
  void evalModel( 
    const Vector            &x
    ,bool                   newx
    ,const ZeroOrderInfo    *zero_order_info  // != NULL if only zero-order info
    ,const ObjGradInfo      *obj_grad_info    // != NULL if obj-grad and below info
    ) const;

};	// end class NLPThyraModelEvaluatorBase

}	// end namespace NLPInterfacePack

#endif	// NLPIP_NLP_THYRA_MODEL_EVALUATOR_BASE_HPP
