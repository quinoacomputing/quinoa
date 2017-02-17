//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_AdjointModelEvaluator.hpp"
#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_Types.hpp"

#include "Thyra_DetachedVectorView.hpp"

#include "../SinCos/SinCosModel.hpp"
#include "../VanderPol/VanderPolModel.hpp"

namespace Rythmos  {

TEUCHOS_UNIT_TEST( Rythmos_AdjointModelEvaluator, create ) {
  RCP<AdjointModelEvaluator<double> > ame = rcp(new AdjointModelEvaluator<double>());
  TEST_ASSERT(!is_null(ame));
}

TEUCHOS_UNIT_TEST( Rythmos_AdjointModelEvaluator, evalLinearModel ) {
  RCP<SinCosModel> fwdModel = sinCosModel(true);
  TimeRange<double> fwdTimeRange(0.0,1.0);
  RCP<AdjointModelEvaluator<double> > ame = adjointModelEvaluator<double>(fwdModel,fwdTimeRange);
  //ame->setVerbLevel(Teuchos::VERB_EXTREME);

  // in args:
  double lambda_0 = 3.0;
  double lambda_1 = 4.0;
  double lambda_dot_0 = 5.0;
  double lambda_dot_1 = 6.0;
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = ame->createInArgs();
  {
    RCP<const Thyra::VectorSpaceBase<double> > space = ame->get_x_space();
    RCP<VectorBase<double> > x = Thyra::createMember(space);
    RCP<VectorBase<double> > xdot = Thyra::createMember(space);
    {
      Thyra::DetachedVectorView<double> x_view( *x );
      x_view[0] = lambda_0;
      x_view[1] = lambda_1;
      Thyra::DetachedVectorView<double> xdot_view( *xdot );
      xdot_view[0] = lambda_dot_0;
      xdot_view[1] = lambda_dot_1;
    }
    inArgs.set_x(x);
    inArgs.set_x_dot(xdot);
  }
  double alpha = 7.0;
  double beta = 8.0;
  inArgs.set_alpha(alpha);
  inArgs.set_beta(beta);
  double time = 0.1;
  inArgs.set_t(time); // doesn't do anything.

  // out args:
  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = ame->createOutArgs();
  RCP<VectorBase<double> > f_out;
  RCP<Thyra::LinearOpBase<double> > W_op;
  {
    RCP<const Thyra::VectorSpaceBase<double> > space = ame->get_f_space();
    f_out = Thyra::createMember(space);
    V_S(Teuchos::outArg(*f_out),0.0);
    outArgs.set_f(f_out);
    W_op = ame->create_W_op(); 
    outArgs.set_W_op(W_op);
  }

  // eval model:
  ame->evalModel(inArgs,outArgs);

  // verify output:
  {
    double f = 1.0;
    double L = 1.0;
    Thyra::ConstDetachedVectorView<double> f_view( *f_out );
    TEST_EQUALITY( f_view[0], lambda_dot_0 + pow(f/L,2.0)*lambda_1 ); // 5.0 + 1.0*4.0 = 9.0
    TEST_EQUALITY( f_view[1], lambda_dot_1 - lambda_0 ); // 6.0 - 3.0 = 3.0
    RCP<const Thyra::VectorSpaceBase<double> > space = ame->get_x_space();
    RCP<VectorBase<double> > x0 = Thyra::createMember(space);
    RCP<VectorBase<double> > x1 = Thyra::createMember(space);
    {
      Thyra::DetachedVectorView<double> x0_view( *x0 );
      x0_view[0] = 1.0;
      x0_view[1] = 0.0;
    }
    Thyra::apply<double>( *W_op, Thyra::NOTRANS, *x0, Teuchos::outArg(*x1) );
    {
      Thyra::ConstDetachedVectorView<double> x1_view( *x1 );
      TEST_EQUALITY( x1_view[0], alpha );
      TEST_EQUALITY( x1_view[1], -beta );
    }
    {
      Thyra::DetachedVectorView<double> x0_view( *x0 );
      x0_view[0] = 0.0;
      x0_view[1] = 1.0;
    }
    Thyra::apply<double>( *W_op, Thyra::NOTRANS, *x0, Teuchos::outArg(*x1) );
    {
      Thyra::ConstDetachedVectorView<double> x1_view( *x1 );
      TEST_EQUALITY( x1_view[0], beta*pow(f/L,2.0) );
      TEST_EQUALITY( x1_view[1], alpha );
    }
  }
}

/*
TEUCHOS_UNIT_TEST( Rythmos_AdjointModelEvaluator, evalNonLinearModel ) {
  RCP<VanderPolModel> fwdModel = vanderPolModel(true);
  TimeRange<double> fwdTimeRange(0.0,1.0);
  RCP<AdjointModelEvaluator<double> > ame = adjointModelEvaluator<double>(fwdModel,fwdTimeRange);
  //ame->setVerbLevel(Teuchos::VERB_EXTREME);

  // in args:
  double lambda_0 = 3.0;
  double lambda_1 = 4.0;
  double lambda_dot_0 = 5.0;
  double lambda_dot_1 = 6.0;
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = ame->createInArgs();
  {
    RCP<const Thyra::VectorSpaceBase<double> > space = ame->get_x_space();
    RCP<VectorBase<double> > x = Thyra::createMember(space);
    RCP<VectorBase<double> > xdot = Thyra::createMember(space);
    {
      Thyra::DetachedVectorView<double> x_view( *x );
      x_view[0] = lambda_0;
      x_view[1] = lambda_1;
      Thyra::DetachedVectorView<double> xdot_view( *xdot );
      xdot_view[0] = lambda_dot_0;
      xdot_view[1] = lambda_dot_1;
    }
    inArgs.set_x(x);
    inArgs.set_x_dot(xdot);
  }
  double alpha = 7.0;
  double beta = 8.0;
  inArgs.set_alpha(alpha);
  inArgs.set_beta(beta);
  double time = 0.1;
  inArgs.set_t(time); // doesn't do anything.

  // out args:
  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = ame->createOutArgs();
  RCP<VectorBase<double> > f_out;
  RCP<Thyra::LinearOpBase<double> > W_op;
  {
    RCP<const Thyra::VectorSpaceBase<double> > space = ame->get_f_space();
    f_out = Thyra::createMember(space);
    V_S(Teuchos::outArg(*f_out),0.0);
    outArgs.set_f(f_out);
    W_op = ame->create_W_op(); 
    outArgs.set_W_op(W_op);
  }

  // eval model:
  ame->evalModel(inArgs,outArgs);

  // verify output:
  {
    double f = 1.0;
    double L = 1.0;
    Thyra::ConstDetachedVectorView<double> f_view( *f_out );
    TEST_EQUALITY( f_view[0], lambda_dot_0 + pow(f/L,2.0)*lambda_1 ); // 5.0 + 1.0*4.0 = 9.0
    TEST_EQUALITY( f_view[1], lambda_dot_1 - lambda_0 ); // 6.0 - 3.0 = 3.0
    RCP<const Thyra::VectorSpaceBase<double> > space = ame->get_x_space();
    RCP<VectorBase<double> > x0 = Thyra::createMember(space);
    RCP<VectorBase<double> > x1 = Thyra::createMember(space);
    {
      Thyra::DetachedVectorView<double> x0_view( *x0 );
      x0_view[0] = 1.0;
      x0_view[1] = 0.0;
    }
    Thyra::apply<double>( *W_op, Thyra::NOTRANS, *x0, Teuchos::outArg(*x1) );
    {
      Thyra::ConstDetachedVectorView<double> x1_view( *x1 );
      TEST_EQUALITY( x1_view[0], alpha );
      TEST_EQUALITY( x1_view[1], -beta );
    }
    {
      Thyra::DetachedVectorView<double> x0_view( *x0 );
      x0_view[0] = 0.0;
      x0_view[1] = 1.0;
    }
    Thyra::apply<double>( *W_op, Thyra::NOTRANS, *x0, Teuchos::outArg(*x1) );
    {
      Thyra::ConstDetachedVectorView<double> x1_view( *x1 );
      TEST_EQUALITY( x1_view[0], beta*pow(f/L,2.0) );
      TEST_EQUALITY( x1_view[1], alpha );
    }
  }
}
*/

} // namespace Rythmos



