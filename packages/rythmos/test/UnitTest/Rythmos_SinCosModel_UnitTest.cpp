//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
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

#include "../SinCos/SinCosModel.hpp"

#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Rythmos {

typedef ModelEvaluatorBase MEB;
using Teuchos::as;

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, setImplicitFlag ) {
  {
    RCP<SinCosModel> explicit_model = rcp(new SinCosModel);
    MEB::InArgs<double> explicit_model_ic;
    TEST_THROW( explicit_model_ic = explicit_model->getNominalValues(), std::logic_error );
    explicit_model->setImplicitFlag(false);
    explicit_model_ic = explicit_model->getNominalValues();
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_beta), true );
  }
  {
    RCP<SinCosModel> implicit_model = rcp(new SinCosModel);
    MEB::InArgs<double> implicit_model_ic;
    TEST_THROW( implicit_model_ic = implicit_model->getNominalValues(), std::logic_error );
    implicit_model->setImplicitFlag(true);
    implicit_model_ic = implicit_model->getNominalValues();
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_beta), true );
  }
}


TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, nominalValues ) {
  double tol = 1.0e-10;
  {
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::InArgs<double> explicit_ic = explicit_model->getNominalValues();
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_beta), true );

    TEST_EQUALITY_CONST( explicit_ic.get_t(), 0.0 );
    RCP<const VectorBase<double> > explicit_ic_x = explicit_ic.get_x();
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*explicit_ic_x,0)+1.0, 1.0, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*explicit_ic_x,1), 1.0, tol );
    TEST_EQUALITY_CONST( explicit_ic.get_beta(), 0.0 );
  }

  {
    RCP<SinCosModel> explicit_model = sinCosModel();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Implicit model formulation",false);
    pl->set("Accept model parameters",true);
    pl->set("Coeff a",25.0);
    pl->set("Coeff f",10.0);
    pl->set("Coeff L",3.0);
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> explicit_ic = explicit_model->getNominalValues();
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_beta), true );

    TEST_EQUALITY_CONST( explicit_ic.get_t(), 0.0 );

    RCP<const VectorBase<double> > explicit_ic_x = explicit_ic.get_x();
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*explicit_ic_x, 0)+1.0, 1.0, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*explicit_ic_x, 1), 1.0, tol );

    TEST_EQUALITY_CONST( explicit_ic.get_beta(), 0.0 );

    RCP<const VectorBase<double> > explicit_ic_p = explicit_ic.get_p(0);
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_p,0), 25.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_p,1), 10.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_p,2),  3.0 );
  }
  
  {
    RCP<SinCosModel> model = sinCosModel();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Implicit model formulation",false);
    pl->set("Provide nominal values",false);
    model->setParameterList(pl);
    MEB::InArgs<double> model_ic = model->getNominalValues();
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_beta), true );

    TEST_EQUALITY_CONST( model_ic.get_t(), 0.0 );
    RCP<const VectorBase<double> > x = model_ic.get_x();
    TEST_EQUALITY_CONST( is_null(x), true );
    TEST_EQUALITY_CONST( model_ic.get_beta(), 0.0 );
  }

  {
    RCP<SinCosModel> implicit_model = sinCosModel(true);
    MEB::InArgs<double> implicit_ic = implicit_model->getNominalValues();
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_t), true);
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_beta), true );

    TEST_EQUALITY_CONST( implicit_ic.get_t(), 0.0 );
    RCP<const VectorBase<double> > implicit_ic_x = implicit_ic.get_x();
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*implicit_ic_x,0), 0.0, tol );
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x,1), 1.0 );
    RCP<const VectorBase<double> > implicit_ic_x_dot = implicit_ic.get_x_dot();
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x_dot,0), 1.0 );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*implicit_ic_x_dot,1), 0.0, tol );
    TEST_EQUALITY_CONST( implicit_ic.get_alpha(), 0.0 );
    TEST_EQUALITY_CONST( implicit_ic.get_beta(), 0.0 );
  }

  {
    RCP<SinCosModel> model = sinCosModel();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Implicit model formulation", true);
    pl->set("Provide nominal values",false);
    model->setParameterList(pl);

    MEB::InArgs<double> model_ic = model->getNominalValues();
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_t), true);
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( model_ic.supports(MEB::IN_ARG_beta), true );

    TEST_EQUALITY_CONST( model_ic.get_t(), 0.0 );
    RCP<const VectorBase<double> > x = model_ic.get_x();
    TEST_EQUALITY_CONST( is_null(x), true );
    RCP<const VectorBase<double> > x_dot = model_ic.get_x_dot();
    TEST_EQUALITY_CONST( is_null(x_dot), true );
    TEST_EQUALITY_CONST( model_ic.get_alpha(), 0.0 );
    TEST_EQUALITY_CONST( model_ic.get_beta(), 0.0 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, p_names ) {
  RCP<SinCosModel> model = sinCosModel();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Accept model parameters", true);
  model->setParameterList(pl);
  RCP<const Teuchos::Array<std::string> > p_names;
#ifdef HAVE_RYTHMOS_DEBUG
  TEST_THROW( p_names = model->get_p_names(1), std::logic_error );
#endif // HAVE_RYTHMOS_DEBUG
  p_names = model->get_p_names(0);
  TEST_EQUALITY_CONST( Teuchos::as<int>(p_names->size()), 3 );
  TEST_EQUALITY_CONST( (*p_names)[0], "Model Coefficient:  a" );
  TEST_EQUALITY_CONST( (*p_names)[1], "Model Coefficient:  f" );
  TEST_EQUALITY_CONST( (*p_names)[2], "Model Coefficient:  L" );
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, spaces ) {
  {
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    RCP<const Thyra::VectorSpaceBase<double> > x_space = explicit_model->get_x_space();
    TEST_EQUALITY_CONST( x_space->dim(), 2 );
    RCP<const Thyra::VectorSpaceBase<double> > f_space = explicit_model->get_f_space();
    TEST_EQUALITY_CONST( f_space->dim(), 2 );
    RCP<const Thyra::VectorSpaceBase<double> > g_space = explicit_model->get_g_space(0);
    TEST_EQUALITY_CONST( g_space->dim(), 1 );
#ifdef HAVE_RYTHMOS_DEBUG
    TEST_THROW( explicit_model->get_g_space(1), std::logic_error );
#endif // HAVE_RYTHMOS_DEBUG
    RCP<const Thyra::VectorSpaceBase<double> > p_space = explicit_model->get_p_space(0);
    TEST_EQUALITY_CONST( is_null(p_space), true );
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    explicit_model->setParameterList(pl);
    p_space = explicit_model->get_p_space(0);
    TEST_EQUALITY_CONST( p_space->dim(), 3 );
#ifdef HAVE_RYTHMOS_DEBUG
    TEST_THROW( explicit_model->get_p_space(1), std::logic_error );
#endif // HAVE_RYTHMOS_DEBUG
  }
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Implicit model formulation",true);
    RCP<SinCosModel> implicit_model = sinCosModel(true);
    RCP<const Thyra::VectorSpaceBase<double> > x_space = implicit_model->get_x_space();
    TEST_EQUALITY_CONST( x_space->dim(), 2 );
    RCP<const Thyra::VectorSpaceBase<double> > f_space = implicit_model->get_f_space();
    TEST_EQUALITY_CONST( f_space->dim(), 2 );
    RCP<const Thyra::VectorSpaceBase<double> > g_space = implicit_model->get_g_space(0);
    TEST_EQUALITY_CONST( g_space->dim(), 1 );
#ifdef HAVE_RYTHMOS_DEBUG
    TEST_THROW( implicit_model->get_g_space(1), std::logic_error );
#endif // HAVE_RYTHMOS_DEBUG
    RCP<const Thyra::VectorSpaceBase<double> > p_space = implicit_model->get_p_space(0);
    TEST_EQUALITY_CONST( is_null(p_space), true );
    pl->set("Accept model parameters",true);
    implicit_model->setParameterList(pl);
    p_space = implicit_model->get_p_space(0);
    TEST_EQUALITY_CONST( p_space->dim(), 3 );
#ifdef HAVE_RYTHMOS_DEBUG
    TEST_THROW( implicit_model->get_p_space(1), std::logic_error );
#endif // HAVE_RYTHMOS_DEBUG
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, create_W_op ) {
  RCP<SinCosModel> explicit_model = sinCosModel(false);
  RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
  RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(matrix), false );
  TEST_EQUALITY_CONST( matrix->domain()->dim(), 2 );
  TEST_EQUALITY_CONST( matrix->range()->dim(), 2 );
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, get_W_factory ) {
  RCP<SinCosModel> explicit_model = sinCosModel(false);
  RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > W_factory = explicit_model->get_W_factory();
  RCP<const Thyra::DefaultSerialDenseLinearOpWithSolveFactory<double> > myFactory =
    Teuchos::rcp_dynamic_cast<const Thyra::DefaultSerialDenseLinearOpWithSolveFactory<double> >(W_factory,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(myFactory), false );
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, createInArgs ) {
  {
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), true );
  }
  {
    RCP<SinCosModel> implicit_model = sinCosModel(true);
    MEB::InArgs<double> inArgs = implicit_model->createInArgs();
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), true );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, createOutArgs ) {
  {
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), true );
  }
  {
    RCP<SinCosModel> explicit_model = sinCosModel(true);
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), true );
  }

}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, exactSolution ) {
  std::vector<double> t_values;
  int N = 10;
  double t = 25;
  for (int i=0 ; i<2*N+1 ; ++i) {
    t_values.push_back( -t + 2*t*i/(2*N) );
  }

  double tol = 1.0e-10;

  double a = 1.0; // [0.0]
  double f = 3.0; // [1.0]
  double L = 4.0; // [1.0]
  double x0 = 0.5; // [0.0]
  double x1 = 1.4; // [1.0]
  double t0 = 0.333; // [0.0]
  double phi = atan(((f/L)/x1)*(x0-a))-(f/L)*t0; 
  double b = x1/((f/L)*cos((f/L)*t0+phi));
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Coeff a", a);
  pl->set("Coeff f", f);
  pl->set("Coeff L", L);
  pl->set("IC x_0", x0);
  pl->set("IC x_1", x1);
  pl->set("IC t_0", t0);
  {
    RCP<SinCosModel> explicit_model = sinCosModel();
    pl->set("Implicit model formulation", false);
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> exact_sol = explicit_model->getExactSolution(0.0);
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_beta), true );
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      MEB::InArgs<double> exact_sol2 = explicit_model->getExactSolution(t_values[i]);
      TEST_FLOATING_EQUALITY( exact_sol2.get_t(), t_values[i], tol );
      RCP<const VectorBase<double> > x = exact_sol2.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,0), a+b*sin((f/L)*t_values[i]+phi), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,1), b*(f/L)*cos((f/L)*t_values[i]+phi), tol );
      TEST_EQUALITY_CONST( exact_sol2.get_beta(), 0.0 );
    }
  }
  {
    RCP<SinCosModel> implicit_model = sinCosModel();
    pl->set("Implicit model formulation", true);
    implicit_model->setParameterList(pl);
    MEB::InArgs<double> exact_sol2 = implicit_model->getExactSolution(0.0);
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_beta), true );
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      MEB::InArgs<double> exact_sol3 = implicit_model->getExactSolution(t_values[i]);
      TEST_FLOATING_EQUALITY( exact_sol3.get_t(), t_values[i], tol );
      RCP<const VectorBase<double> > x = exact_sol3.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,0), a+b*sin((f/L)*t_values[i]+phi), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,1), b*(f/L)*cos((f/L)*t_values[i]+phi), tol );
      RCP<const VectorBase<double> > x_dot = exact_sol3.get_x_dot();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x_dot,0), b*(f/L)*cos((f/L)*t_values[i]+phi), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x_dot,1), -b*(f/L)*(f/L)*sin((f/L)*t_values[i]+phi), tol );
      TEST_EQUALITY_CONST( exact_sol3.get_alpha(), 0.0 );
      TEST_EQUALITY_CONST( exact_sol3.get_beta(), 0.0 );
    }
  }
}

/* x_0(t) = a + b*sin((f/L)*t+phi)
 * x_1(t) = b*(f/L)*cos((f/L)*t+phi)
 *
 * xdot_0(t) = b*(f/L)*cos((f/L)*t+phi)
 * xdot_1(t) = -b*(f/L)^2*sin((f/L)*t+phi)
 *
 * p = (a,f,L)
 */
TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, exactSensSolution ) {
  std::vector<double> t_values;
  int N = 10;
  {
    double t = 25;
    for (int i=0 ; i<2*N+1 ; ++i) {
      t_values.push_back( -t + 2*t*i/(2*N) );
    }
  }

  double tol = 1.0e-10;

  double a = 1.0; // [0.0]
  double f = 3.0; // [1.0]
  double L = 4.0; // [1.0]
  double x0 = 0.5; // [0.0]
  double x1 = 1.4; // [1.0]
  double t0 = 0.333; // [0.0]
  double phi = atan(((f/L)/x1)*(x0-a))-(f/L)*t0; 
  double b = x1/((f/L)*cos((f/L)*t0+phi));
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Accept model parameters",true);
  pl->set("Coeff a", a);
  pl->set("Coeff f", f);
  pl->set("Coeff L", L);
  pl->set("IC x_0", x0);
  pl->set("IC x_1", x1);
  pl->set("IC t_0", t0);
  {
    RCP<SinCosModel> explicit_model = sinCosModel();
    pl->set("Implicit model formulation", false);
    explicit_model->setParameterList(pl);
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      double t = t_values[i];
      MEB::InArgs<double> exact_sol = explicit_model->getExactSensSolution(0,t);
      RCP<const VectorBase<double> > sens = exact_sol.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,0), 1.0, tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,1), 0.0, tol );

      exact_sol = explicit_model->getExactSensSolution(1,t_values[i]);
      sens = exact_sol.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,0), (b*t/L)*cos((f/L)*t+phi), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,1), (b/L)*cos((f/L)*t+phi)-(b*f*t/(L*L))*sin((f/L)*t+phi), tol );

      exact_sol = explicit_model->getExactSensSolution(2,t_values[i]);
      sens = exact_sol.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,0), -(b*f*t/(L*L))*cos((f/L)*t+phi), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,1), -(b*f/(L*L))*cos((f/L)*t+phi)+(b*f*f*t/(L*L*L))*sin((f/L)*t+phi), tol );
    }
  }
  {
    RCP<SinCosModel> implicit_model = sinCosModel();
    pl->set("Implicit model formulation", true);
    implicit_model->setParameterList(pl);
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      double t = t_values[i];
      MEB::InArgs<double> exact_sol = implicit_model->getExactSensSolution(0,t);
      RCP<const VectorBase<double> > sens = exact_sol.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,0), 1.0, tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,1), 0.0, tol );
      RCP<const VectorBase<double> > sens_dot = exact_sol.get_x_dot();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens_dot,0), 0.0, tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens_dot,1), 0.0, tol );

      exact_sol = implicit_model->getExactSensSolution(1,t);
      sens = exact_sol.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,0), (b*t/L)*cos((f/L)*t+phi), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,1), (b/L)*cos((f/L)*t+phi)-(b*f*t/(L*L))*sin((f/L)*t+phi), tol );
      sens_dot = exact_sol.get_x_dot();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens_dot,0), (b/L)*cos((f/L)*t+phi)-(b*f*t/(L*L))*sin((f/L)*t+phi),  tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens_dot,1), -(2.0*b*f/(L*L))*sin((f/L)*t+phi)-(b*f*f*t/(L*L*L))*cos((f/L)*t+phi), tol );

      exact_sol = implicit_model->getExactSensSolution(2,t);
      sens = exact_sol.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,0), -(b*f*t/(L*L))*cos((f/L)*t+phi), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,1), -(b*f/(L*L))*cos((f/L)*t+phi)+(b*f*f*t/(L*L*L))*sin((f/L)*t+phi), tol );
      sens_dot = exact_sol.get_x_dot();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens_dot,0), -(b*f/(L*L))*cos((f/L)*t+phi)+(b*f*f*t/(L*L*L))*sin((f/L)*t+phi),  tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens_dot,1), +(2.0*b*f*f/(L*L*L))*sin((f/L)*t+phi)+(b*f*f*f*t/(L*L*L*L))*cos((f/L)*t+phi), tol );
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, evalExplicitModel ) {
  double a = 2.0;
  double freq = 3.0;
  double L = 4.0;
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Implicit model formulation", false);
  pl->set("Coeff a", a);
  pl->set("Coeff f", freq);
  pl->set("Coeff L", L);
  { // Explicit model, just load f.
    RCP<SinCosModel> explicit_model = sinCosModel();
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    double t = 4.0; 
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);

    RCP<VectorBase<double> > f = Thyra::createMember(explicit_model->get_f_space());
    outArgs.set_f(f);

    explicit_model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), (freq/L)*(freq/L)*(a-5.0) );
  }
  { // Explicit model, load f and W_op
    RCP<SinCosModel> explicit_model = sinCosModel();
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    double t = 4.0; 
    double beta = 7.0;
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_beta(beta);

    RCP<VectorBase<double> > f = Thyra::createMember(explicit_model->get_f_space());
    RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
    outArgs.set_f(f);
    outArgs.set_W_op(W_op);

    explicit_model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), (freq/L)*(freq/L)*(a-5.0) );

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 0.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), 7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,0), -(freq/L)*(freq/L)*7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,1), 0.0 );
  }
  { // Explicit model, just load W_op
    RCP<SinCosModel> explicit_model = sinCosModel();
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    double t = 4.0; 
    double beta = 7.0;
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_beta(beta);

    RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
    outArgs.set_W_op(W_op);

    explicit_model->evalModel(inArgs,outArgs);

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 0.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), 7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,0), -(freq/L)*(freq/L)*7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,1), 0.0 );
  }

}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, evalImplicitModel ) {
  double a = 2.0;
  double freq = 3.0;
  double L = 4.0;
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Implicit model formulation", true);
  pl->set("Coeff a", a);
  pl->set("Coeff f", freq);
  pl->set("Coeff L", L);
  { // Implicit model, just load f.
    RCP<SinCosModel> implicit_model = sinCosModel();
    implicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = implicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = implicit_model->createOutArgs();
    double t = 4.0; 
    RCP<VectorBase<double> > x = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = 8.0;
      x_dot_view[1] = 9.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);

    RCP<VectorBase<double> > f = Thyra::createMember(implicit_model->get_f_space());
    outArgs.set_f(f);

    implicit_model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 8.0 - 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), 9.0 - (freq/L)*(freq/L)*(a-5.0) );
  }
  { // Implicit model, load f and W_op
    RCP<SinCosModel> implicit_model = sinCosModel();
    implicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = implicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = implicit_model->createOutArgs();
    double t = 4.0; 
    double alpha = 11.0;
    double beta = 7.0;
    RCP<VectorBase<double> > x = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = 8.0;
      x_dot_view[1] = 9.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);
    inArgs.set_alpha(alpha);
    inArgs.set_beta(beta);

    RCP<VectorBase<double> > f = Thyra::createMember(implicit_model->get_f_space());
    RCP<Thyra::LinearOpBase<double> > W_op = implicit_model->create_W_op();
    outArgs.set_f(f);
    outArgs.set_W_op(W_op);

    implicit_model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 8.0 - 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), 9.0 - (freq/L)*(freq/L)*(a-5.0) );

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 11.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), -7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,0), (freq/L)*(freq/L)*7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,1), 11.0 );
  }
  { // Implicit model, just load W_op
    RCP<SinCosModel> implicit_model = sinCosModel();
    implicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = implicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = implicit_model->createOutArgs();
    double t = 4.0; 
    double alpha = 11.0;
    double beta = 7.0;
    RCP<VectorBase<double> > x = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = 8.0;
      x_dot_view[1] = 9.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);
    inArgs.set_alpha(alpha);
    inArgs.set_beta(beta);

    RCP<Thyra::LinearOpBase<double> > W_op = implicit_model->create_W_op();
    outArgs.set_W_op(W_op);

    implicit_model->evalModel(inArgs,outArgs);

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 11.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), -7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,0), (freq/L)*(freq/L)*7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,1), 11.0 );
  }

}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, modelParams ) {
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Accept model parameters", true);
  { // Explicit model, just load f.
    RCP<SinCosModel> model = sinCosModel();
    pl->set("Implicit model formulation", false);
    model->setParameterList(pl);
    MEB::InArgs<double> inArgs = model->createInArgs();
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    double t = 4.0; 
    RCP<VectorBase<double> > x = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    double a = 25.0;
    double freq = 10.0;
    double L = 3.0;
    RCP<VectorBase<double> > p = Thyra::createMember(model->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p_view(*p);
      p_view[0] = a;
      p_view[1] = freq;
      p_view[2] = L;
    }
    inArgs.set_p(0,p);

    RCP<VectorBase<double> > f = Thyra::createMember(model->get_f_space());
    outArgs.set_f(f);

    model->evalModel(inArgs,outArgs);

    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), 6.0, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), (freq/L)*(freq/L)*(a-5.0), tol );
  }
  { // Implicit model, just load f.
    RCP<SinCosModel> model = sinCosModel();
    pl->set("Implicit model formulation", true);
    model->setParameterList(pl);
    MEB::InArgs<double> inArgs = model->createInArgs();
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    double t = 4.0; 
    RCP<VectorBase<double> > x = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    RCP<VectorBase<double> > xdot = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> xdot_view(*xdot);
      xdot_view[0] = 8.0;
      xdot_view[1] = 9.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(xdot);
    double a = 25.0;
    double freq = 10.0;
    double L = 3.0;
    RCP<VectorBase<double> > p = Thyra::createMember(model->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p_view(*p);
      p_view[0] = a;
      p_view[1] = freq;
      p_view[2] = L;
    }
    inArgs.set_p(0,p);

    RCP<VectorBase<double> > f = Thyra::createMember(model->get_f_space());
    outArgs.set_f(f);

    model->evalModel(inArgs,outArgs);

    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), 8.0-6.0, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), 9.0-(freq/L)*(freq/L)*(a-5.0), tol );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, DfDp ) {
  double a = 1.0;
  double freq = 2.0;
  double L = 3.0;
  double x0 = 1.2;
  double x1 = 1.3;
  //double t0 = 0.0;
  //double phi = atan(((freq/L)/x1)*(x0-a))-(freq/L)*t0; 
  //double b = x1/((freq/L)*cos((freq/L)*t0+phi));

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Accept model parameters", true);
  pl->set("Coeff a", a);
  pl->set("Coeff f", freq);
  pl->set("Coeff L", L);
  pl->set("IC x_0", x0);
  pl->set("IC x_1", x1);
  { // Explicit model, load f and DfDp
    RCP<SinCosModel> model = sinCosModel();
    pl->set("Implicit model formulation", false);
    model->setParameterList(pl);
    MEB::InArgs<double> inArgs = model->createInArgs();
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    double t = 4.0; 
    RCP<VectorBase<double> > x = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    double a2 = 25.0;
    double freq2 = 10.0;
    double L2 = 4.0;
    RCP<VectorBase<double> > p = Thyra::createMember(model->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p_view(*p);
      p_view[0] = a2;
      p_view[1] = freq2;
      p_view[2] = L2;
    }
    inArgs.set_p(0,p);

    RCP<VectorBase<double> > f = Thyra::createMember(model->get_f_space());
    outArgs.set_f(f);

    RCP<Thyra::MultiVectorBase<double> > mv = Thyra::createMembers(model->get_f_space(),3);
    MEB::Derivative<double> DfDp(mv);
    outArgs.set_DfDp(0, DfDp);

    model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), (freq2/L2)*(freq2/L2)*(a2-5.0) );

    // Check that DfDp is correct
    {
      Thyra::ConstDetachedMultiVectorView<double> mv_view( *mv );
      TEST_EQUALITY_CONST( mv_view(0,0), 0.0 ); // df0dp0 = df0da = 0.0
      TEST_EQUALITY_CONST( mv_view(0,1), 0.0 ); // df0dp1 = df0df = 0.0
      TEST_EQUALITY_CONST( mv_view(0,2), 0.0 ); // df0dp2 = df0dL = 0.0
      TEST_EQUALITY_CONST( mv_view(1,0), (freq2/L2)*(freq2/L2) ); // df1dp0 = df1da = (f/L2)*(f/L2)
      TEST_EQUALITY_CONST( mv_view(1,1), (2.0*freq2/(L2*L2))*(a2-5.0) ); // df1dp1 = df1df = 2*f/(L2*L2)*(a2-x0)
      TEST_EQUALITY_CONST( mv_view(1,2), -(2.0*freq2*freq2/(L2*L2*L2))*(a2-5.0) ); // df1dp2 = df1dL = -2*f*f/(L2*L2*L2)*(a2-x0)
    }
    
  }
  { // Implicit model, load f and DfDp
    RCP<SinCosModel> model = sinCosModel();
    pl->set("Implicit model formulation", true);
    model->setParameterList(pl);
    MEB::InArgs<double> inArgs = model->createInArgs();
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    double t = 4.0; 
    RCP<VectorBase<double> > x = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = 8.0;
      x_dot_view[1] = 9.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);
    double a2 = 25.0;
    double freq2 = 10.0;
    double L2 = 4.0;
    RCP<VectorBase<double> > p = Thyra::createMember(model->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p_view(*p);
      p_view[0] = a2;
      p_view[1] = freq2;
      p_view[2] = L2;
    }
    inArgs.set_p(0,p);

    RCP<VectorBase<double> > f = Thyra::createMember(model->get_f_space());
    outArgs.set_f(f);

    RCP<Thyra::MultiVectorBase<double> > mv = Thyra::createMembers(model->get_f_space(),3);
    MEB::Derivative<double> DfDp(mv);
    outArgs.set_DfDp(0, DfDp);

    model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 8.0 - 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), 9.0 - (freq2/L2)*(freq2/L2)*(a2-5.0) );

    // Check that DfDp is correct
    {
      Thyra::ConstDetachedMultiVectorView<double> mv_view( *mv );
      TEST_EQUALITY_CONST( mv_view(0,0), 0.0 );
      TEST_EQUALITY_CONST( mv_view(0,1), 0.0 );
      TEST_EQUALITY_CONST( mv_view(0,2), 0.0 );
      TEST_EQUALITY_CONST( mv_view(1,0), -(freq2/L2)*(freq2/L2) );
      TEST_EQUALITY_CONST( mv_view(1,1), -(2.0*freq2/(L2*L2))*(a2-5.0) );
      TEST_EQUALITY_CONST( mv_view(1,2), +(2.0*freq2*freq2/(L2*L2*L2))*(a2-5.0) );
    }
    
  }
}

} // namespace Rythmos



