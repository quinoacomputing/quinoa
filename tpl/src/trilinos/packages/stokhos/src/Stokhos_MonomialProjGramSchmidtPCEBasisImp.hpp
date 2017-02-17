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

#include "Stokhos_SDMUtils.hpp"
#include "Stokhos_OrthogonalizationFactory.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
MonomialProjGramSchmidtPCEBasis(
  ordinal_type max_p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::ParameterList& params) :
  GSReducedPCEBasisBase<ordinal_type,value_type>(max_p, pce, quad, params),
  name("Monomial Proj Gram Schmidt PCE Basis")
{
  this->setup(max_p, pce, quad);
}

template <typename ordinal_type, typename value_type>
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
~MonomialProjGramSchmidtPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
buildReducedBasis(
  ordinal_type max_p, 
  value_type threshold,
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& A, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& F,
  const Teuchos::Array<value_type>& weights, 
  Teuchos::Array< Stokhos::MultiIndex<ordinal_type> >& terms_,
  Teuchos::Array<ordinal_type>& num_terms_,
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& Qp_, 
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& Q_)
{
  // Compute basis terms -- 2-D array giving powers for each linear index
  ordinal_type max_sz;
  CPBUtils::compute_terms(max_p, this->d, max_sz, terms_, num_terms_);

  // Compute B matrix -- monomials in F
  // for i=0,...,nqp-1
  //   for j=0,...,sz-1
  //      B(i,j) = F(i,1)^terms_[j][1] * ... * F(i,d)^terms_[j][d]
  // where sz is the total size of a basis up to order p and terms_[j] 
  // is an array of powers for each term in the total-order basis
  ordinal_type nqp = weights.size();
  SDM B(nqp, max_sz);
  for (ordinal_type i=0; i<nqp; i++) {
    for (ordinal_type j=0; j<max_sz; j++) {
      B(i,j) = 1.0;
      for (ordinal_type k=0; k<this->d; k++)
	B(i,j) *= std::pow(F(i,k), terms_[j][k]);
    }
  }

  // Project B into original basis -- should use SPAM for this
  SDM Bp(this->pce_sz, max_sz);
  const Teuchos::Array<value_type>& basis_norms = 
    this->pce_basis->norm_squared();
  for (ordinal_type i=0; i<this->pce_sz; i++) {
    for (ordinal_type j=0; j<max_sz; j++) {
      Bp(i,j) = 0.0;
      for (ordinal_type k=0; k<nqp; k++)
	Bp(i,j) += weights[k]*B(k,j)*A(k,i);
      Bp(i,j) /= basis_norms[i];
    }
  }

  // Rescale columns of Bp to have unit norm
  for (ordinal_type j=0; j<max_sz; j++) {
    value_type nrm = 0.0;
    for (ordinal_type i=0; i<this->pce_sz; i++)
      nrm += Bp(i,j)*Bp(i,j)*basis_norms[i];
    nrm = std::sqrt(nrm);
    for (ordinal_type i=0; i<this->pce_sz; i++)
      Bp(i,j) /= nrm;
  }

  // Compute our new basis -- each column of Qp is the coefficients of the
  // new basis in the original basis.  Constraint pivoting so first d+1
  // columns and included in Qp.
  Teuchos::Array<value_type> w(this->pce_sz, 1.0);
  SDM R;
  Teuchos::Array<ordinal_type> piv(max_sz);
  for (int i=0; i<this->d+1; i++)
    piv[i] = 1;
  typedef Stokhos::OrthogonalizationFactory<ordinal_type,value_type> SOF;
  ordinal_type sz_ = SOF::createOrthogonalBasis(
    this->orthogonalization_method, threshold, this->verbose, Bp, w, 
    Qp_, R, piv);

  // Evaluate new basis at original quadrature points
  Q_.reshape(nqp, sz_);
  ordinal_type ret = 
    Q_.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, Qp_, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  return sz_;
}
