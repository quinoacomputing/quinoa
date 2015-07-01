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

#ifndef NLPIP_NLP_FIRST_ORDER_THYRA_MODEL_EVALUATOR_HPP
#define NLPIP_NLP_FIRST_ORDER_THYRA_MODEL_EVALUATOR_HPP

#include <vector>

#include "NLPInterfacePack_NLPThyraModelEvaluatorBase.hpp"

namespace NLPInterfacePack {

/** \brief Implement the %NLPFirstOrder interface using a
 * <tt>Thyra::ModelEvaluator</tt> object.
 *
 * ToDo: Finish documentation!
 */
class NLPFirstOrderThyraModelEvaluator
  : virtual public NLPFirstOrder
  , virtual public NLPThyraModelEvaluatorBase
{
public:

  /** \brief Initialize to uninitialized */
  NLPFirstOrderThyraModelEvaluator();

  /** \brief Calls <tt>initialize()</tt>. */
  NLPFirstOrderThyraModelEvaluator(
    const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
    );

  /** \brief .Initialize given a <tt>Thyra::ModelEvaluator</tt> and
   * a description of how to interpret it.
   *
   * ToDo: Finish documentation!
   *
   * Todo: Add arguments for auxiliary inequalites and equalities
   */
  void initialize(
    const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
    );

  /** @name Overridden public members from NLP */
  //@{

  /** \brief . */
  void initialize(bool test_setup);
  /** \brief . */
  void unset_quantities();

  //@}

  /** @name Overridden public members from NLPFirstOrder */
  //@{

  /** \brief Overridden to check the concrete type of Gc */
  void set_Gc(MatrixOp* Gc);
  /** \brief . */
  const NLPFirstOrder::mat_fcty_ptr_t factory_Gc() const;
  /** \brief Returns an ExampleBasisSystem */
  const basis_sys_ptr_t basis_sys() const;

  //@}

protected:

  /** @name Overridden protected members from NLPFirstOrder */
  //@{

  /** \brief . */
  void imp_calc_Gc(
    const Vector& x, bool newx
    ,const FirstOrderInfo& first_order_info) const;

  //@}

private:

  // /////////////////////////////////////////
  // Private member functions

  /** \brief . */
  void evalModel( 
    const Vector            &x
    ,bool                   newx
    ,const ZeroOrderInfo    *zero_order_info  // != NULL if only zero-order info
    ,const ObjGradInfo      *obj_grad_info    // != NULL if obj-grad and below info
    ,const FirstOrderInfo   *first_order_info // != NULL if first-order and below info
    ) const;

};	// end class NLPFirstOrderThyraModelEvaluator

}	// end namespace NLPInterfacePack

#endif	// NLPIP_NLP_FIRST_ORDER_THYRA_MODEL_EVALUATOR_HPP
