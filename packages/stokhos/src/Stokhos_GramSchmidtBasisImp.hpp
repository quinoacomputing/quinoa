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

#include "Teuchos_BLAS.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
GramSchmidtBasis(
   const Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& basis_,
   const Teuchos::Array< Teuchos::Array<value_type> >& points,
   const Teuchos::Array<value_type>& weights_,
   const value_type& sparse_tol_) :
  name("Gram Schmidt Basis"),
  basis(basis_),
  weights(weights_),
  sparse_tol(sparse_tol_),
  p(basis->order()),
  d(basis->dimension()),
  sz(basis->size()),
  norms(sz),
  gs_mat(sz,sz),
  basis_vals_tmp(sz)
{
  // Get quadrature data
  ordinal_type nqp = weights.size();
  Teuchos::Array< Teuchos::Array<value_type> > values(nqp);
  for (ordinal_type k=0; k<nqp; k++) {
    values[k].resize(sz);
    basis->evaluateBases(points[k], values[k]);
  }

  // Compute all inner products
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> inner_product(sz,sz);
  inner_product.putScalar(0.0);
  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<=i; j++) {
      value_type t = 0.0;
      for (ordinal_type k=0; k<nqp; k++)
	t += weights[k]*values[k][i]*values[k][j];
      inner_product(i,j) = t;
    }
  }

  // Classical Gram-Schmidt algorithm:
  // u_i = v_i - \sum_{j<i} (v_i,u_j)/(u_j,u_j) u_j
  // u_j = \sum_{k<=i} a_{jk} v_k
  // => u_i = v_i - \sum_{j<i}\sum_{k<=j} (v_i,u_j)/(u_j,u_j)*a_{jk}*v_k
  for (ordinal_type i=0; i<sz; i++) {

    // a_{ii} = 1.0
    gs_mat(i,i) = 1.0;

    for (ordinal_type j=0; j<i; j++) {

      // compute t = (v_i,u_j)/(u_j,u_j)
      value_type t = 0.0;
      for (ordinal_type k=0; k<=j; k++)
	t += gs_mat(j,k)*inner_product(i,k);
      t /= norms[j];
      
      // substract contribution to a_{ik}:  t*a_{jk}
      for (ordinal_type k=0; k<=j; k++)
	gs_mat(i,k) -= t*gs_mat(j,k);
    }

    // compute (u_i,u_i) = \sum_{j,k<=i} a_{ij}*a_{ik}*(v_j,v_k)
    value_type nrm = 0.0;
    for (ordinal_type j=0; j<=i; j++) {
      for (ordinal_type k=0; k<=j; k++)
	nrm += gs_mat(i,j)*gs_mat(i,k)*inner_product(j,k);
      for (ordinal_type k=j+1; k<=i; k++)
      	nrm += gs_mat(i,j)*gs_mat(i,k)*inner_product(k,j);
    }
    norms[i] = nrm;

  }

  basis_values.resize(nqp);
  for (ordinal_type k=0; k<nqp; k++) {
    basis_values[k].resize(sz);
    for (ordinal_type i=0; i<sz; i++) {
      value_type t = 0.0;
      for (ordinal_type j=0; j<=i; j++)
	t += gs_mat(i,j)*values[k][j];
      basis_values[k][i] = t;
    }
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
~GramSchmidtBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
computeTripleProductTensor() const
{
  Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
    Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>);
  ordinal_type nqp = weights.size();
  for (ordinal_type j=0; j<sz; j++) {
    for (ordinal_type i=0; i<sz; i++) {
      for (ordinal_type k=0; k<sz; k++) {
	value_type t = 0.0;
	for (ordinal_type l=0; l<nqp; l++)
	  t += 
	    weights[l]*basis_values[l][i]*basis_values[l][j]*basis_values[l][k];
	if (std::abs(t) > sparse_tol)
	  Cijk->add_term(i,j,k,t);
      }
    }
  }

  Cijk->fillComplete();

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
computeLinearTripleProductTensor() const
{
  Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
    Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>);
  ordinal_type nqp = weights.size();
  for (ordinal_type j=0; j<sz; j++) {
    for (ordinal_type i=0; i<sz; i++) {
      for (ordinal_type k=0; k<d+1; k++) {
	value_type t = 0.0;
	for (ordinal_type l=0; l<nqp; l++)
	  t += 
	    weights[l]*basis_values[l][i]*basis_values[l][j]*basis_values[l][k];
	if (std::abs(t) > sparse_tol)
	  Cijk->add_term(i,j,k,t);
      }
    }
  }

  Cijk->fillComplete();

  return Cijk;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  value_type z = 0.0;
  for (ordinal_type j=0; j<sz; j++)
    z += gs_mat(i,j)*basis->evaluateZero(j);

  return z;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
evaluateBases(const Teuchos::ArrayView<const value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  basis->evaluateBases(point, basis_vals_tmp);
  for (ordinal_type i=0; i<sz; i++) {
    value_type t = 0.0;
    for (ordinal_type j=0; j<sz; j++)
      t += gs_mat(i,j)*basis_vals_tmp[j];
    basis_vals[i] = t;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Gram-Schmidt basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Matrix coefficients:\n";
  os << gs_mat << std::endl;
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<sz; i++)
    os << norms[i] << " ";
  os << "\n";
  os << "Underlying basis:\n";
  os << *basis;
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
transformCoeffs(const value_type *in, value_type *out) const
{
  Teuchos::BLAS<ordinal_type, value_type> blas;
  for (ordinal_type i=0; i<sz; i++)
    out[i] = in[i];
  blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::TRANS,
	    Teuchos::UNIT_DIAG, sz, 1, 1.0, gs_mat.values(), sz, out, sz);
}
