/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#ifndef THYRA_SIMPLE_2D_MODEL_EVALUATOR_DEF_HPP
#define THYRA_SIMPLE_2D_MODEL_EVALUATOR_DEF_HPP


#include "Thyra_Simple2DModelEvaluator_decl.hpp"
#include "Thyra_SimpleDenseLinearOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Thyra {


// Nonmember constuctors


template<class Scalar>
Teuchos::RCP<Simple2DModelEvaluator<Scalar> >
simple2DModelEvaluator()
{
  return Teuchos::rcp(new Simple2DModelEvaluator<Scalar>);
}


// Initializers/Accessors


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::set_d(const Scalar &d)
{
  d_ = d;
}


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::set_p(const Teuchos::ArrayView<const Scalar> &p)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(p_.size(), p.size());
#endif
  p_().assign(p);
}


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<Scalar> x0(x0_);
  x0.sv().values()().assign(x0_in);
}


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::setShowGetInvalidArgs(bool showGetInvalidArg)
{
  showGetInvalidArg_ = showGetInvalidArg;
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Simple2DModelEvaluator<Scalar>::get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Simple2DModelEvaluator<Scalar>::get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Simple2DModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Simple2DModelEvaluator<Scalar>::create_W_op() const
{
  return createNonconstSimpleDenseLinearOp<Scalar>(
    createMembers<Scalar>(f_space_, x_space_->dim())
    );
}


template<class Scalar>
Teuchos::RCP<Thyra::PreconditionerBase<Scalar> >
Simple2DModelEvaluator<Scalar>::create_W_prec() const
{
  return nonconstUnspecifiedPrec<Scalar>(
    createNonconstSimpleDenseLinearOp<Scalar>(
      createMembers<Scalar>(f_space_, x_space_->dim())
      )
    );
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
Simple2DModelEvaluator<Scalar>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Simple2DModelEvaluator<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Simple2DModelEvaluator<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::rcp_dynamic_cast;
  const Scalar one = 1.0, two = 2.0, zero = 0.0;

  const ConstDetachedVectorView<Scalar> x(inArgs.get_x());

  const RCP<Thyra::VectorBase<Scalar> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase< Scalar > > W_op_out = outArgs.get_W_op();
  const RCP<Thyra::PreconditionerBase< Scalar > > W_prec_out = outArgs.get_W_prec();

  if (nonnull(f_out)) {
    const DetachedVectorView<Scalar> f(f_out);
    f[0] = x[0] + x[1] * x[1] - p_[0];
    f[1] = d_ * (x[0] * x[0] - x[1] - p_[1]);
  }

  if (nonnull(W_op_out)) {
    const RCP<SimpleDenseLinearOp<Scalar> > W =
      rcp_dynamic_cast<SimpleDenseLinearOp<Scalar> >(W_op_out, true);
    const RCP<MultiVectorBase<Scalar> > W_mv = W->getNonconstMultiVector();
    Thyra::DetachedMultiVectorView<Scalar> W_dmvv(W_mv);
    W_dmvv(0, 0) = one;
    W_dmvv(0, 1) = two * x[1];
    W_dmvv(1, 0) = d_ * two * x[0];
    W_dmvv(1, 1) = -d_;
  }

  if (nonnull(W_prec_out)) {
    const RCP<SimpleDenseLinearOp<Scalar> > W_prec_op =
      rcp_dynamic_cast<SimpleDenseLinearOp<Scalar> >(
        W_prec_out->getNonconstUnspecifiedPrecOp(), true);
    const RCP<MultiVectorBase<Scalar> > W_prec_mv = W_prec_op->getNonconstMultiVector();
    Thyra::DetachedMultiVectorView<Scalar> W_prec_dmvv(W_prec_mv);
    // Diagonal inverse of W (see W above)
    W_prec_dmvv(0, 0) = one;
    W_prec_dmvv(0, 1) = zero;
    W_prec_dmvv(1, 0) = zero;
    W_prec_dmvv(1, 1) = -one/d_;
  }
  
}


// private


template<class Scalar>
Simple2DModelEvaluator<Scalar>::Simple2DModelEvaluator()
  : x_space_(Thyra::defaultSpmdVectorSpace<Scalar>(2)),
    f_space_(x_space_),
    W_factory_(Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>()),
    d_(0.0),
    p_(x_space_->dim(), Teuchos::ScalarTraits<Scalar>::zero()),
    showGetInvalidArg_(false)
{

  using Teuchos::RCP;
  using Thyra::VectorBase;
  using Thyra::createMember;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  prototypeInArgs_ = inArgs;
  
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.setSupports(MEB::OUT_ARG_W_prec);
  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  x0_ = createMember(x_space_);
  V_S(x0_.ptr(), ST::zero());
  nominalValues_.set_x(x0_);

  set_d(10.0);
  set_p(Teuchos::tuple<Scalar>(2.0, 0.0)());
  set_x0(Teuchos::tuple<Scalar>(1.0, 1.0)());

}


} // namespace Thyra


//
// Explicit instantiation macro
//
// Must be expanded from within the global namespace!
//

#define SIMPLE_2D_MODEL_EVALUATOR_INSTANT(SCALAR) \
  \
  template class Simple2DModelEvaluator<SCALAR >; \
  \
  template Teuchos::RCP<Simple2DModelEvaluator<SCALAR > > \
  simple2DModelEvaluator(); \


#endif // THYRA_SIMPLE_2D_MODEL_EVALUATOR_DEF_HPP
