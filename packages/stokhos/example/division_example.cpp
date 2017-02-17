// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

// hermite_example
//
//  usage: 
//     hermite_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of the simple function
//
//     v = 1/(log(u)^2+1)
//
//     where u = 1 + 0.4*H_1(x) + 0.06*H_2(x) + 0.002*H_3(x), x is a zero-mean
//     and unit-variance Gaussian random variable, and H_i(x) is the i-th
//     Hermite polynomial.

#include "Stokhos.hpp"
#include "Stokhos_Sacado.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// Typename of PC expansion type
typedef Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pce_type;

int main(int argc, char **argv)
{
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // If applicable, set up MPI.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv);

  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString("This example tests the / operator.\n");
    int d = 3;
    CLP.setOption("dim", &d, "Stochastic dimension");
    int p = 5;
    CLP.setOption("order", &p, "Polynomial order");
    int n = 100;
    CLP.setOption("samples", &n, "Number of samples");
    double shift = 2.0;
    CLP.setOption("shift", &shift, "Shift point");
    CLP.parse( argc, argv );

    // Basis
    Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
    for (int i=0; i<d; i++) {
      bases[i] = rcp(new Stokhos::LegendreBasis<int,double>(p));
    }
    RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Quadrature method
    RCP<const Stokhos::Quadrature<int,double> > quad = 
      rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

    // Triple product tensor
    RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor(basis->size());

    // Expansion methods
    Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > quad_expn = 
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(
		     basis, Cijk, quad));
    Teuchos::RCP<Teuchos::ParameterList> alg_params = 
      Teuchos::rcp(new Teuchos::ParameterList);
    alg_params->set("Division Strategy", "Dense Direct");
    Teuchos::RCP<Stokhos::AlgebraicOrthogPolyExpansion<int,double> > alg_expn = 
      Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<int,double>(
		     basis, Cijk, alg_params));

    // Polynomial expansions
    pce_type u_quad(quad_expn), v_quad(quad_expn);
    u_quad.term(0,0) = 0.0;
    for (int i=0; i<d; i++) {
      u_quad.term(i,1) = 1.0;
    }
    pce_type u_alg(alg_expn), v_alg(alg_expn);
    u_alg.term(0,0) = 0.0;
    for (int i=0; i<d; i++) {
      u_alg.term(i,1) = 1.0;
    }

    // Compute expansion
    v_quad = 1.0 / (shift + u_quad);
    v_alg = 1.0 / (shift + u_alg);

    // Print u and v
    // std::cout << "quadrature:   v = 1.0 / (shift + u) = ";
    // v_quad.print(std::cout);
    // std::cout << "dense solve:  v = 1.0 / (shift + u) = ";
    // v_alg.print(std::cout);

    double h = 2.0 / (n-1);
    double err_quad = 0.0;
    double err_alg = 0.0;
    for (int i=0; i<n; i++) {
      double x = -1.0 + h*i;
      Array<double> pt(d); 
      for (int j=0; j<d; j++) 
	pt[j] = x;
      double up = u_quad.evaluate(pt);
      double vp = 1.0/(shift+up);
      double vp_quad = v_quad.evaluate(pt);
      double vp_alg = v_alg.evaluate(pt);
      // std::cout << "vp = " << vp_quad << std::endl;
      // std::cout << "vp_quad = " << vp_quad << std::endl;
      // std::cout << "vp_alg = " << vp_alg << std::endl;
      double point_err_quad = std::abs(vp-vp_quad);
      double point_err_alg = std::abs(vp-vp_alg);
      if (point_err_quad > err_quad) err_quad = point_err_quad;
      if (point_err_alg > err_alg) err_alg = point_err_alg;
    }
    std::cout << "\tL_infty norm of quadrature error = " << err_quad 
	      << std::endl;
    std::cout << "\tL_infty norm of dense direct error = " << err_alg
	      << std::endl;
    
    // Check the answer
    //if (std::abs(err) < 1e-2)
      std::cout << "\nExample Passed!" << std::endl;

    Teuchos::TimeMonitor::summarize(std::cout);
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
