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

#include "Stokhos_ReducedQuadratureFactory.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
GSReducedPCEBasisBase(
  ordinal_type max_p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::ParameterList& params_) :
  params(params_),
  pce_basis(pce[0].basis()),
  pce_sz(pce_basis->size()),
  p(max_p),
  d(pce.size()),
  verbose(params.get("Verbose", false)),
  rank_threshold(params.get("Rank Threshold", 1.0e-12)),
  orthogonalization_method(params.get("Orthogonalization Method", 
				      "Householder"))
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
setup(
  ordinal_type max_p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad)
{
  // Check for pce's that are constant and don't represent true random
  // dimensions
  Teuchos::Array< const Stokhos::OrthogPolyApprox<ordinal_type, value_type>* > pce2;
  for (ordinal_type i=0; i<pce.size(); i++) {
    if (pce[i].standard_deviation() > 1.0e-15)
      pce2.push_back(&pce[i]);
  }
  d = pce2.size();

  // Get quadrature data
  const Teuchos::Array<value_type>& weights = quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<value_type> >& points = 
    quad->getQuadPoints(); 
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values = 
    quad->getBasisAtQuadPoints();
  ordinal_type nqp = weights.size();

  // Original basis at quadrature points -- needed to transform expansions
  // in this basis back to original
  SDM A(nqp, pce_sz);
  for (ordinal_type i=0; i<nqp; i++)
    for (ordinal_type j=0; j<pce_sz; j++)
      A(i,j) = basis_values[i][j];

  // Compute norms of each pce for rescaling
  Teuchos::Array<value_type> pce_norms(d, 0.0);
  for (ordinal_type j=0; j<d; j++) {
    for (ordinal_type i=0; i<pce_sz; i++)
      pce_norms[j] += (*pce2[j])[i]*(*pce2[j])[i]*pce_basis->norm_squared(i);
    pce_norms[j] = std::sqrt(pce_norms[j]);
  }

  // Compute F matrix -- PCEs evaluated at all quadrature points
  // Since F is used in the reduced quadrature below as the quadrature points
  // for this reduced basis, does scaling by the pce_norms mess up the points?
  // No -- F essentially defines the random variables this basis is a function
  // of, and thus they can be scaled in any way we want.  Because we don't 
  // explicitly write the basis in terms of F, the scaling is implicit.
  SDM F(nqp, d);
  Teuchos::Array< Teuchos::Array<value_type> > values(nqp);
  for (ordinal_type i=0; i<nqp; i++) 
    for (ordinal_type j=0; j<d; j++)
      F(i,j) = pce2[j]->evaluate(points[i], basis_values[i]);

  // Build the reduced basis
  sz = buildReducedBasis(max_p, rank_threshold, A, F, weights, terms, num_terms,
			 Qp, Q);

  // Compute reduced quadrature rule
  Teuchos::ParameterList quad_params = params.sublist("Reduced Quadrature");
  Stokhos::ReducedQuadratureFactory<ordinal_type,value_type> quad_factory(
    quad_params);
  SDM Q2;
  if (quad_params.isParameter("Reduced Quadrature Method") &&
      quad_params.get<std::string>("Reduced Quadrature Method") == "Q2") {
    Teuchos::Array< Stokhos::MultiIndex<ordinal_type> > terms2;
    Teuchos::Array<ordinal_type> num_terms2;
    value_type rank_threshold2 = quad_params.get("Q2 Rank Threshold", 
						 rank_threshold);
    SDM Qp2;
    //ordinal_type sz2 = 
    buildReducedBasis(2*max_p, rank_threshold2, A, F, weights, terms2, 
		      num_terms2, Qp2, Q2);
  }
  reduced_quad = quad_factory.createReducedQuadrature(Q, Q2, F, weights);

  // Basis is orthonormal by construction
  norms.resize(sz, 1.0);
}

template <typename ordinal_type, typename value_type>
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
~GSReducedPCEBasisBase()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
computeTripleProductTensor() const

{
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
computeLinearTripleProductTensor() const

{
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
evaluateBases(const Teuchos::ArrayView<const value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Gram-Schmidt basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Matrix coefficients:\n";
  os << Qp << std::endl;
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<sz; i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
transformToOriginalBasis(const value_type *in, value_type *out,
			 ordinal_type ncol, bool transpose) const
{
  if (transpose) {
    SDM zbar(Teuchos::View, const_cast<value_type*>(in), ncol, ncol, sz);
    SDM z(Teuchos::View, out, ncol, ncol, pce_sz);
    ordinal_type ret = 
      z.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, zbar, Qp, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
  else {
    SDM zbar(Teuchos::View, const_cast<value_type*>(in), sz, sz, ncol);
    SDM z(Teuchos::View, out, pce_sz, pce_sz, ncol);
    ordinal_type ret = 
      z.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Qp, zbar, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
transformFromOriginalBasis(const value_type *in, value_type *out,
			 ordinal_type ncol, bool transpose) const
{
  if (transpose) {
    SDM z(Teuchos::View, const_cast<value_type*>(in), ncol, ncol, pce_sz);
    SDM zbar(Teuchos::View, out, ncol, ncol, sz);
    ordinal_type ret = 
      zbar.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, z, Qp, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
  else {
    SDM z(Teuchos::View, const_cast<value_type*>(in), pce_sz, pce_sz, ncol);
    SDM zbar(Teuchos::View, out, sz, sz, ncol);
    ordinal_type ret = 
      zbar.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Qp, z, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
Stokhos::GSReducedPCEBasisBase<ordinal_type, value_type>::
getReducedQuadrature() const
{
  return reduced_quad;
}
