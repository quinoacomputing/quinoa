// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2009) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"
#include "Stokhos_MonomialGramSchmidtSimplexPCEBasis.hpp"

typedef Stokhos::LegendreBasis<int,double> basis_type;

int main(int argc, char **argv)
{
  try {
 
    const unsigned int d = 2;
    const unsigned int p = 5;

    // Create product basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (unsigned int i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new basis_type(p, true));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    std::cout << "original basis size = " << basis->size() << std::endl;

    // Create approximation
    Stokhos::OrthogPolyApprox<int,double> x(basis), u(basis), v(basis), 
      w(basis), w2(basis);
    for (unsigned int i=0; i<d; i++) {
      x.term(i, 1) = 1.0;
    }

    // Tensor product quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    // Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
    //   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis,
    // 								 p+1));

    std::cout << "original quadrature size = " << quad->size() << std::endl;

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor(basis->size());
    
    // Quadrature expansion
    Stokhos::QuadOrthogPolyExpansion<int,double> quad_exp(basis, Cijk, quad);
    
    // Compute PCE via quadrature expansion
    quad_exp.sin(u,x);
    quad_exp.exp(v,x);
    quad_exp.times(w,v,u);
    
    // Create new basis from u and v
    Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > pces(2);
    pces[0] = u;
    pces[1] = v;
    Teuchos::ParameterList params;
    //params.set("Verbose", true);
    //params.set("Reduced Quadrature Method", "L1 Minimization");
    params.set("Reduced Quadrature Method", "Column-Pivoted QR");
    //params.set("Orthogonalization Method", "Classical Gram-Schmidt");
    params.set("Orthogonalization Method", "Modified Gram-Schmidt");
    Teuchos::RCP< Stokhos::MonomialGramSchmidtSimplexPCEBasis<int,double> > gs_basis = 
      Teuchos::rcp(new Stokhos::MonomialGramSchmidtSimplexPCEBasis<int,double>(
		     p, pces, quad, params));
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > gs_quad =
      gs_basis->getReducedQuadrature();
    Stokhos::OrthogPolyApprox<int,double>  u_gs(gs_basis), v_gs(gs_basis), 
      w_gs(gs_basis);
    gs_basis->computeTransformedPCE(0, u_gs);
    gs_basis->computeTransformedPCE(1, v_gs);

    std::cout << "reduced basis size = " << gs_basis->size() << std::endl;
    std::cout << "reduced quadrature size = " << gs_quad->size() << std::endl;

    Stokhos::OrthogPolyApprox<int,double> u2(basis), v2(basis); 
    gs_basis->transformToOriginalBasis(u_gs.coeff(), u2.coeff());
    gs_basis->transformToOriginalBasis(v_gs.coeff(), v2.coeff());
    double err_u = 0.0;
    double err_v = 0.0;
    for (int i=0; i<basis->size(); i++) {
      double eu = std::abs(u[i]-u2[i]);
      double ev = std::abs(v[i]-v2[i]);
      if (eu > err_u) err_u = eu;
      if (ev > err_v) err_v = ev;
    }
    std::cout << "error in u transformation = " << err_u << std::endl;
    std::cout << "error in v transformation = " << err_v << std::endl;
    
    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > gs_Cijk =
      Teuchos::null;
    
    // Gram-Schmidt quadrature expansion
    Teuchos::RCP< Teuchos::ParameterList > gs_exp_params = 
      Teuchos::rcp(new Teuchos::ParameterList);
    gs_exp_params->set("Use Quadrature for Times", true);
    Stokhos::QuadOrthogPolyExpansion<int,double> gs_quad_exp(gs_basis, 
							     gs_Cijk,
							     gs_quad,
							     gs_exp_params);
    
    // Compute w_gs = u_gs*v_gs in Gram-Schmidt basis
    gs_quad_exp.times(w_gs, u_gs, v_gs);
    
    // Project w_gs back to original basis
    gs_basis->transformToOriginalBasis(w_gs.coeff(), w2.coeff());

    std::cout.precision(12);
    // std::cout << "w = " << std::endl << w;
    // std::cout << "w2 = " << std::endl << w2;
    // std::cout << "w_gs = " << std::endl << w_gs;

    double err_w = 0.0;
    for (int i=0; i<basis->size(); i++) {
      double ew = std::abs(w[i]-w2[i]);
      if (ew > err_w) err_w = ew;
    }
    
    std::cout.setf(std::ios::scientific);
    std::cout << "w.mean()       = " << w.mean() << std::endl
	      << "w2.mean()      = " << w2.mean() << std::endl
	      << "mean error     = " 
	      << std::abs(w.mean()-w2.mean()) << std::endl
	      << "w.std_dev()    = " << w.standard_deviation() << std::endl
	      << "w2.std_dev()   = " << w2.standard_deviation() << std::endl
	      << "std_dev error  = " 
	      << std::abs(w.standard_deviation()-w2.standard_deviation()) 
	      << std::endl
	      << "w coeff error  = " << err_w << std::endl;
    
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
