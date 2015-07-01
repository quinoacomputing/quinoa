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

#ifndef NLPIP_NLP_DIRECT_THYRA_MODEL_EVALUATOR_HPP
#define NLPIP_NLP_DIRECT_THYRA_MODEL_EVALUATOR_HPP

#include "NLPInterfacePack_NLPThyraModelEvaluatorBase.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace NLPInterfacePack {

/** \brief Implement the %NLPFirstOrder interface using a
 * <tt>Thyra::ModelEvaluator</tt> object.
 *
 * ToDo: Finish documentation!
 */
class NLPDirectThyraModelEvaluator
  : virtual public NLPDirect
  , virtual public NLPThyraModelEvaluatorBase
{
public:

  /** \brief Utility object that computes directional finite differences for objective */
  STANDARD_COMPOSITION_MEMBERS( Thyra::DirectionalFiniteDiffCalculator<value_type>, objDirecFiniteDiffCalculator );

  /** \brief Utility object that computes directional finite differences for constraints */
  STANDARD_COMPOSITION_MEMBERS( Thyra::DirectionalFiniteDiffCalculator<value_type>, conDirecFiniteDiffCalculator );

  /** \brief Set if model.DfDp is constant or not */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, DfDp_is_const );

  /** \brief Initialize to uninitialized */
  NLPDirectThyraModelEvaluator();

  /** \brief Calls <tt>initialize()</tt>. */
  NLPDirectThyraModelEvaluator(
    const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
    ,const objDirecFiniteDiffCalculator_ptr_t                       objDirecFiniteDiffCalculator = Teuchos::null
    ,const conDirecFiniteDiffCalculator_ptr_t                       conDirecFiniteDiffCalculator = Teuchos::null
    );

  /** \brief .Initialize given a <tt>Thyra::ModelEvaluator</tt> and
   * a description of how to interpret it.
   *
   * ToDo: Finish documentation!
   */
  void initialize(
    const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
    ,const objDirecFiniteDiffCalculator_ptr_t                       objDirecFiniteDiffCalculator = Teuchos::null
    ,const conDirecFiniteDiffCalculator_ptr_t                       conDirecFiniteDiffCalculator = Teuchos::null
    );

  /** @name Overridden public members from NLP */
  //@{

  /** \brief . */
  void initialize(bool test_setup);
  /** \brief . */
  void unset_quantities();

  //@}

  /** @name Overridden public members from NLPObjGrad */
  //@{

  /** \brief . */
  bool supports_Gf() const;
  /** \brief . */
  bool supports_Gf_prod() const;
  /** \brief . */
  value_type calc_Gf_prod(const Vector& x, const Vector& d, bool newx) const;

  //@}

  /** @name Overridden public members from NLPDirect */
  //@{

  /** \brief . */
  Range1D var_dep() const;
  /** \brief . */
  Range1D var_indep() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_D() const;
  /** \brief . */
  const mat_sym_nonsing_fcty_ptr_t factory_S() const;
  /** \brief . */
  void calc_point(
    const Vector     &x
    ,value_type      *f
    ,VectorMutable   *c
    ,bool            recalc_c
    ,VectorMutable   *Gf
    ,VectorMutable   *py
    ,VectorMutable   *rGf
    ,MatrixOp        *GcU
    ,MatrixOp        *D
    ,MatrixOp        *Uz
    ) const;
  /** \brief . */
  void calc_semi_newton_step(
    const Vector    &x
    ,VectorMutable  *c
    ,bool           recalc_c
    ,VectorMutable  *py
    ) const;

  //@}

private:

  mutable Teuchos::RCP<Thyra::LinearOpWithSolveBase<value_type> >  thyra_C_;
  mutable Teuchos::RCP<Thyra::MultiVectorBase<value_type> >        thyra_N_;
  
};	// end class NLPDirectThyraModelEvaluator

}	// end namespace NLPInterfacePack

#endif	// NLPIP_NLP_DIRECT_THYRA_MODEL_EVALUATOR_HPP
