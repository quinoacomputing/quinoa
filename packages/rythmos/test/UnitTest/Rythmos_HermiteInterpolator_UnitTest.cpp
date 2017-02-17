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

#include "Rythmos_HermiteInterpolator.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

#include "Teuchos_Polynomial.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

using Teuchos::Polynomial;
using Teuchos::RCP;

TEUCHOS_UNIT_TEST( Rythmos_HermiteInterpolator, nonMemberConstructor ) {
  RCP<HermiteInterpolator<double> > hi = hermiteInterpolator<double>();
  TEST_EQUALITY_CONST( hi->order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_HermiteInterpolator, interpolate ) {
  RCP<InterpolatorBase<double> > hi = hermiteInterpolator<double>();
  double maxOrder = hi->order();
  for (int order = 0 ; order <= maxOrder+1 ; ++order) {
    TEST_EQUALITY( order, order );
    Polynomial<double> p(order,1.0); 
    RCP<DataStore<double>::DataStoreVector_t> data_in = rcp( new DataStore<double>::DataStoreVector_t );
    double T0 = 0.0;
    double T1 = 1.0;
    int N = 5;
    for (int i=0 ; i < N ; ++i) {
      double t = ((T1-T0)/(N-1.0))*i+T0;
      double x = 0.0;
      double xdot = 0.0;
      p.evaluate(t,&x,&xdot);
      RCP<Thyra::VectorBase<double> > xv, xvdot;
      double accuracy = 0.0;
      xv = createDefaultVector<double>(1,x);
      xvdot = createDefaultVector<double>(1,xdot);
      data_in->push_back(DataStore<double>(t,xv,xvdot,accuracy));
    }
    Array<double> t_values;
    DataStore<double>::DataStoreVector_t data_out;
    N = 2*N;
    for (int i=0 ; i < N ; ++i) {
      double t = ((T1-T0)/(N-1.0))*i+T0;
      t_values.push_back(t);
    }
    interpolate<double>(*hi, data_in, t_values, &data_out);
    // Verify that the interpolated values are exactly the same as the polynomial values
    // unless the order of polynomial is greater than the order of the interpolator
    unsigned int N_out = data_out.size();
    for (unsigned int i=0 ; i < N_out ; ++i ) {
      double x = 0.0;
      double xdot = 0.0;
      double t = data_out[i].time;
      RCP<const Thyra::VectorBase<double> > xv = data_out[i].x;
      RCP<const Thyra::VectorBase<double> > xvdot = data_out[i].xdot;
      {
        Thyra::ConstDetachedVectorView<double> xv_view(*xv);
        x = xv_view[0];
        Thyra::ConstDetachedVectorView<double> xvdot_view(*xvdot);
        xdot = xvdot_view[0];
      }
      double x_exact = 0.0;
      double xdot_exact = 0.0;
      p.evaluate(t,&x_exact,&xdot_exact);
      double tol = 1.0e-15;
      if ((order <= maxOrder) || (i == 0) || (i == N_out-1)) {
        TEST_FLOATING_EQUALITY( x, x_exact, tol );
        TEST_FLOATING_EQUALITY( xdot, xdot_exact, tol );
      } else {
        TEST_COMPARE( fabs((x-x_exact)/x_exact), >, tol );
        TEST_COMPARE( fabs((xdot-xdot_exact)/xdot_exact), >, tol );
      }
    }
  }
}


} // namespace Rythmos



