// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Stokhos_Sacado.hpp"
#include "Stokhos_Sacado_Kokkos.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// The function to compute the polynomial chaos expansion of,
// written as a template function
template <class ScalarType>
ScalarType simple_function(const ScalarType& u) {
  ScalarType z = std::log(u);
  return 1.0/(z*z + 1.0);
}

int main(int argc, char **argv)
{
  // Typename of Polynomial Chaos scalar type
  typedef Stokhos::StandardStorage<int,double> pce_storage_type;
  typedef Sacado::ETPCE::OrthogPoly<double, pce_storage_type> pce_type;

  // Typename of ensemble scalar type
  const int EnsembleSize = 8;
  typedef Stokhos::StaticFixedStorage<int,double,EnsembleSize,Kokkos::DefaultExecutionSpace> ensemble_storage_type;
  typedef Sacado::MP::Vector<ensemble_storage_type> ensemble_type;

  // Short-hand for several classes used below
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Stokhos::OneDOrthogPolyBasis;
  using Stokhos::HermiteBasis;
  using Stokhos::LegendreBasis;
  using Stokhos::CompletePolynomialBasis;
  using Stokhos::Quadrature;
  using Stokhos::TotalOrderIndexSet;
  using Stokhos::SmolyakSparseGridQuadrature;
  using Stokhos::TensorProductQuadrature;
  using Stokhos::Sparse3Tensor;
  using Stokhos::QuadOrthogPolyExpansion;

  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example computes the PC expansion of a simple function.\n");
    int p = 4;
    CLP.setOption("order", &p, "Polynomial order");
    bool sparse = false;
    CLP.setOption("sparse", "tensor", &sparse,
                  "Use sparse grid or tensor product quadrature");

    // Parse arguments
    CLP.parse( argc, argv );

    // Basis of dimension 3, order given by command-line option
    const int d = 3;
    Array< RCP<const OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++) {
      bases[i] = rcp(new HermiteBasis<int,double>(p, true));
    }
    RCP<const CompletePolynomialBasis<int,double> > basis =
      rcp(new CompletePolynomialBasis<int,double>(bases));
    const int pce_size = basis->size();
    std::cout << "basis size = " << pce_size << std::endl;

    // Quadrature method
    RCP<const Quadrature<int,double> > quad;
    if (sparse) {
      const TotalOrderIndexSet<int> index_set(d, p);
      quad = rcp(new SmolyakSparseGridQuadrature<int,double>(basis, index_set));
    }
    else {
      quad = rcp(new TensorProductQuadrature<int,double>(basis));
    }
    std::cout << "quadrature size = " << quad->size() << std::endl;

    // Triple product tensor
    RCP<Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();

    // Expansion method
    RCP<QuadOrthogPolyExpansion<int,double> > expn =
      rcp(new QuadOrthogPolyExpansion<int,double>(basis, Cijk, quad));

    // Polynomial expansion of u (note:  these are coefficients in the
    // normalized basis)
    pce_type u(expn);
    u.term(0,0) = 1.0;     // zeroth order term
    u.term(0,1) = 0.1;     // first order term for dimension 0
    u.term(1,1) = 0.05;    // first order term for dimension 1
    u.term(2,1) = 0.01;    // first order term for dimension 2

    //
    // Compute PCE expansion of function using NISP with ensemble propagation
    //

    // Extract quadrature data
    const int num_quad_points                 = quad->size();
    const Array<double>& quad_weights         = quad->getQuadWeights();
    const Array< Array<double> >& quad_points = quad->getQuadPoints();
    const Array< Array<double> >& quad_values = quad->getBasisAtQuadPoints();

    // Loop over quadrature points in blocks of size EnsembleSize
    pce_type v(expn);
    ensemble_type u_ensemble;
    for (int qp_block=0; qp_block<num_quad_points; qp_block+=EnsembleSize) {
      const int qp_sz = qp_block+EnsembleSize <= num_quad_points ?
        EnsembleSize : num_quad_points-qp_block;

      // Evaluate u at each quadrature point
      for (int qp=0; qp<qp_sz; ++qp)
        u_ensemble.fastAccessCoeff(qp) =
          u.evaluate(quad_points[qp_block+qp], quad_values[qp_block+qp]);
      for (int qp=qp_sz; qp<EnsembleSize; ++qp)
        u_ensemble.fastAccessCoeff(qp) = u_ensemble.fastAccessCoeff(qp_sz-1);

      // Evaluate function at each quadrature point
      ensemble_type v_ensemble = simple_function(u_ensemble);

      // Sum results into PCE integral
      for (int pc=0; pc<pce_size; ++pc)
        for (int qp=0; qp<qp_sz; ++qp)
          v.fastAccessCoeff(pc) += v_ensemble.fastAccessCoeff(qp)*quad_weights[qp_block+qp]*quad_values[qp_block+qp][pc];
    }

    /*
    for (int qp=0; qp<num_quad_points; ++qp) {
      double u_qp = u.evaluate(quad_points[qp]);
      double v_qp = simple_function(u_qp);
      double w = quad_weights[qp];
      for (int pc=0; pc<pce_size; ++pc)
        v.fastAccessCoeff(pc) += v_qp*w*quad_values[qp][pc];
    }
    */

    // Print u and v
    std::cout << "\tu = ";
    u.print(std::cout);
    std::cout << "\tv = ";
    v.print(std::cout);

    // Compute moments
    double mean = v.mean();
    double std_dev = v.standard_deviation();

    // Evaluate PCE and function at a point = 0.25 in each dimension
    Teuchos::Array<double> pt(d);
    for (int i=0; i<d; i++)
      pt[i] = 0.25;
    double up = u.evaluate(pt);
    double vp = simple_function(up);
    double vp2 = v.evaluate(pt);

    // Print results
    std::cout << "\tv mean         = " << mean << std::endl;
    std::cout << "\tv std. dev.    = " << std_dev << std::endl;
    std::cout << "\tv(0.25) (true) = " << vp << std::endl;
    std::cout << "\tv(0.25) (pce)  = " << vp2 << std::endl;

     // Check the answer
    if (std::abs(vp - vp2) < 1e-2)
      std::cout << "\nExample Passed!" << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
