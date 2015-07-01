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

#ifndef BARRIER_NLP_H
#define BARRIER_NLP_H

#include "NLPInterfacePack_NLPObjGrad.hpp"

namespace NLPInterfacePack {

/** \brief Simple wrapper that provides an objective fn with the barrier term
 *   included.
 *
 */
class NLPBarrier : public NLPObjGrad
{
public:
    
  /** @name Public Methods */
  //@{

  /// Set the barrier parameter.
  void mu(const value_type mu);

  /// Get the barrier term.  Must be called after <tt>calc_f()</tt>.
  value_type barrier_term() const;

  /// Get the true objective term value.  Must be called after <tt>calc_f()</tt>.
  value_type objective_term() const;

  /// Get the value of the gradient of the barrier term.  Must be called after <tt>calc_Gf()</tt>
  const Teuchos::RCP<Vector> grad_barrier_term() const;

  /// Get the value of the gradient of the true objective term.  Must be called after <tt>calc_Gf()</tt>.
  const Teuchos::RCP<Vector> grad_objective_term() const;
    
  //@}
      
  /** @name Constructors / initializers */
  //@{

  /** \brief Constructor.
   */
  NLPBarrier();

  /** \brief . */
  void InitializeFromNLP(
    Teuchos::RCP<NLP> original_nlp
    );

  //@}

  /** @name Overridden public members from NLPObjGrad */
  //@{

  /** \brief . */
  void initialize(bool test_setup)
  { nlp_->initialize(test_setup); }
  /** \brief . */
  bool is_initialized() const
  { return nlp_->is_initialized(); }
  /** \brief . */
  void set_Gf(VectorMutable* Gf)
  { nlp_->set_Gf(Gf); }
  /** \brief . */
  VectorMutable* get_Gf()
  { return nlp_->get_Gf(); }
  /** \brief . */
  VectorMutable& Gf()
  { return nlp_->Gf(); }
  /** \brief . */
  const Vector& Gf() const
  { return nlp_->Gf(); }
  /// Overloaded to include barrier term
  void calc_Gf(const Vector& x, bool newx = true) const;
  /** \brief . */
  size_type num_Gf_evals() const
  { return nlp_->num_Gf_evals(); }

  //@}

  /** @name Overridden public members from NLP */
  //@{

  /** \brief . */
  void force_xinit_in_bounds(bool force_xinit_in_bounds)
  { nlp_->force_xinit_in_bounds(force_xinit_in_bounds); }
  /** \brief . */
  bool force_xinit_in_bounds() const
  { return nlp_->force_xinit_in_bounds(); }
  /** \brief . */
  size_type n() const
  { return nlp_->n(); }
  /** \brief . */
  size_type m() const
  { return nlp_->m(); }
  /** \brief . */
  vec_space_ptr_t space_x() const
  { return nlp_->space_x(); }
  /** \brief . */
  vec_space_ptr_t space_c() const
  { return nlp_->space_c(); }
  /** \brief . */
  size_type num_bounded_x() const
  { return nlp_->num_bounded_x(); }
  /** \brief . */
  const Vector& xl() const
  { return nlp_->xl(); }
  /** \brief . */
  const Vector& xu() const 
  { return nlp_->xu(); }
  /** \brief . */
  value_type max_var_bounds_viol() const
  { return nlp_->max_var_bounds_viol(); }
  /** \brief . */
  const Vector& xinit() const
  { return nlp_->xinit(); }
  /** \brief . */
  void get_init_lagrange_mult(
    VectorMutable*   lambda
    ,VectorMutable*  nu
    ) const
  { nlp_->get_init_lagrange_mult(lambda, nu); }
  /** \brief . */
  void set_f(value_type* f)
  { nlp_->set_f(f); }
  /** \brief . */
  value_type* get_f()
  { return nlp_->get_f(); }
  /** \brief . */
  value_type& f()
  { return nlp_->f(); }
  /** \brief . */
  const value_type& f() const
  { return nlp_->f(); }
  /** \brief . */
  void set_c(VectorMutable* c)
  { nlp_->set_c(c); }
  /** \brief . */
  VectorMutable* get_c()
  { return nlp_->get_c(); }
  /** \brief . */
  VectorMutable& c()
  { return nlp_->c(); }
  /** \brief . */
  const Vector& c() const
  { return nlp_->c(); }
  /** \brief . */
  void scale_f( value_type scale_f )
  { nlp_->scale_f(); }
  /** \brief . */
  value_type scale_f() const
  { return nlp_->scale_f(); }
  /// Overloaded to include barrier term
  void calc_f(const Vector& x, bool newx =  true) const;
  /** \brief . */
  void calc_c(const Vector& x, bool newx = true) const
  { nlp_->calc_c(x, newx); }
  /** \brief . */
  void report_final_solution(
    const Vector&    x
    ,const Vector*   lambda
    ,const Vector*   nu
    ,bool            is_optimal
    )
  { nlp_->report_final_solution(
    x, lambda, nu, is_optimal
    );
  }
  /** \brief . */
  size_type num_f_evals() const
  { return nlp_->num_f_evals(); }
  /** \brief . */
  size_type num_c_evals() const
  { return nlp_->num_c_evals(); }
  /** \brief . */
  size_type ns() const
  { return nlp_->ns(); }
  /** \brief . */
  vec_space_ptr_t space_c_breve() const
  { return nlp_->space_c_breve(); }
  /** \brief . */
  vec_space_ptr_t space_h_breve() const
  { return nlp_->space_h_breve(); }
  /** \brief . */
  const Vector& hl_breve() const
  { return nlp_->hl_breve(); }
  /** \brief . */
  const Vector& hu_breve() const
  { return nlp_->hu_breve(); }
  /** \brief . */
  void set_c_breve(VectorMutable* c_breve)
  { nlp_->set_c_breve(c_breve); }
  /** \brief . */
  VectorMutable* get_c_breve()
  { return nlp_->get_c_breve(); }
  /** \brief . */
  VectorMutable& c_breve()
  { return nlp_->c_breve(); }
  /** \brief . */
  const Vector& c_breve() const
  { return nlp_->c_breve(); }
  /** \brief . */
  void set_h_breve(VectorMutable* h_breve)
  { nlp_->set_h_breve(h_breve); }
  /** \brief . */
  VectorMutable* get_h_breve()
  { return nlp_->get_h_breve(); }
  /** \brief . */
  VectorMutable& h_breve()
  { return nlp_->h_breve(); }
  /** \brief . */
  const Vector& h_breve() const
  { return nlp_->h_breve(); }
  /** \brief . */
  const Permutation& P_var() const
  { return nlp_->P_var(); }
  /** \brief . */
  const Permutation& P_equ() const
  { return nlp_->P_equ(); }
  /** \brief . */
  void calc_c_breve(const Vector& x, bool newx ) const
  { nlp_->calc_c_breve(x,newx); }
  /** \brief . */
  void calc_h_breve(const Vector& x, bool newx ) const
  { nlp_->calc_h_breve(x,newx); }

  //@}

protected:

  /* protected members Overridden from NLP */
  //@{

  /** \brief . */
  void imp_calc_f(
    const Vector& x
    ,bool newx 
    ,const ZeroOrderInfo& zero_order_info
    ) const;
  /** \brief . */
  void imp_calc_c(
    const Vector& x
    ,bool newx 
    ,const ZeroOrderInfo& zero_order_info
    ) const;
  /** \brief . */
  void imp_calc_c_breve(
    const Vector& x
    ,bool newx 
    ,const ZeroOrderInfo& zero_order_info_breve
    ) const;
  /** \brief . */
  void imp_calc_h_breve(
    const Vector& x
    ,bool newx 
    ,const ZeroOrderInfo& zero_order_info_breve
    ) const;

  //@}

  /* protected members Overridden from NLPObjGrad */
  //@{

  /** \brief . */
  void imp_calc_Gf(
    const Vector& x,
    bool newx, 
    const ObjGradInfo& obj_grad_info
    ) const;

  //@}

private:

  Teuchos::RCP<NLPObjGrad> nlp_;
  value_type                                           mu_;
  mutable value_type                                   barrier_term_;
  mutable value_type                                   objective_term_;
  mutable Teuchos::RCP<VectorMutable>     grad_barrier_term_;
  mutable Teuchos::RCP<VectorMutable>     grad_barrier_term_temp_;
  mutable Teuchos::RCP<VectorMutable>     grad_objective_term_;

  value_type CalculateBarrierTerm(const Vector& x) const;

}; // end class NLPBarrier

} // end namespace NLPInterfacePack

#endif	// BARRIER_NLP_H
