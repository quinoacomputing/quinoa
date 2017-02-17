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
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TestForException.hpp"

template <typename ordinal_type, typename value_type, typename ordering_type>
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
TotalOrderBasis(
  const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > >& bases_,
  const value_type& sparse_tol_,
  const ordering_type& coeff_compare) :
  p(0),
  d(bases_.size()),
  sz(0),
  bases(bases_),
  sparse_tol(sparse_tol_),
  max_orders(d),
  basis_set(coeff_compare),
  norms()
{

  // Compute largest order
  for (ordinal_type i=0; i<d; i++) {
    max_orders[i] = bases[i]->order();
    if (max_orders[i] > p)
      p = max_orders[i];
  }

  // Compute basis terms
  MultiIndex<ordinal_type> orders(d);
  for (ordinal_type i=0; i<d; ++i)
    orders[i] = bases[i]->order();
  AnisotropicTotalOrderIndexSet<ordinal_type> index_set(p, orders);
  ProductBasisUtils::buildProductBasis(index_set, basis_set, basis_map);
  sz = basis_map.size();
    
  // Compute norms
  norms.resize(sz);
  value_type nrm;
  for (ordinal_type k=0; k<sz; k++) {
    nrm = value_type(1.0);
    for (ordinal_type i=0; i<d; i++)
      nrm = nrm * bases[i]->norm_squared(basis_map[k][i]);
    norms[k] = nrm;
  }

  // Create name
  name = "Tensor product basis (";
  for (ordinal_type i=0; i<d-1; i++)
    name += bases[i]->getName() + ", ";
  name += bases[d-1]->getName() + ")";

  // Allocate array for basis evaluation
  basis_eval_tmp.resize(d);
  for (ordinal_type j=0; j<d; j++)
    basis_eval_tmp[j].resize(max_orders[j]+1);
}

template <typename ordinal_type, typename value_type, typename ordering_type>
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
~TotalOrderBasis()
{
}

template <typename ordinal_type, typename value_type, typename ordering_type>
ordinal_type
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
ordinal_type
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
ordinal_type
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
const Teuchos::Array<value_type>&
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
const value_type&
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type, typename ordering_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
computeTripleProductTensor() const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Triple-Product Tensor Fill Time");
#endif

  TotalOrderPredicate<ordinal_type> predicate(p, max_orders);
  
  return ProductBasisUtils::computeTripleProductTensor(
    bases, basis_set, basis_map, predicate, predicate, sparse_tol);
}

template <typename ordinal_type, typename value_type, typename ordering_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
computeLinearTripleProductTensor() const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Triple-Product Tensor Fill Time");
#endif

  TotalOrderPredicate<ordinal_type> predicate(p, max_orders);
  TotalOrderPredicate<ordinal_type> k_predicate(1, max_orders);
  
  return ProductBasisUtils::computeTripleProductTensor(
    bases, basis_set, basis_map, predicate, k_predicate, sparse_tol);
}

template <typename ordinal_type, typename value_type, typename ordering_type>
value_type
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
evaluateZero(ordinal_type i) const
{
  // z = psi_{i_1}(0) * ... * psi_{i_d}(0) where i_1,...,i_d are the basis
  // terms for coefficient i
  value_type z = value_type(1.0);
  for (ordinal_type j=0; j<d; j++)
    z = z * bases[j]->evaluate(value_type(0.0), basis_map[i][j]);

  return z;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
void
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
evaluateBases(const Teuchos::ArrayView<const value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  for (ordinal_type j=0; j<d; j++)
    bases[j]->evaluateBases(point[j], basis_eval_tmp[j]);

  // Only evaluate basis upto number of terms included in basis_pts
  for (ordinal_type i=0; i<sz; i++) {
    value_type t = value_type(1.0);
    for (ordinal_type j=0; j<d; j++)
      t *= basis_eval_tmp[j][basis_map[i][j]];
    basis_vals[i] = t;
  }
}

template <typename ordinal_type, typename value_type, typename ordering_type>
void
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
print(std::ostream& os) const
{
  os << "Tensor product basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Component bases:\n";
  for (ordinal_type i=0; i<d; i++)
    os << *bases[i];
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<static_cast<ordinal_type>(norms.size()); i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type, typename ordering_type>
const Stokhos::MultiIndex<ordinal_type>&
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
term(ordinal_type i) const
{
  return basis_map[i];
}

template <typename ordinal_type, typename value_type, typename ordering_type>
ordinal_type
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
index(const MultiIndex<ordinal_type>& term) const
{
  typename coeff_set_type::const_iterator it = basis_set.find(term);
  TEUCHOS_TEST_FOR_EXCEPTION(it == basis_set.end(), std::logic_error,
			     "Invalid term " << term);
  return it->second;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
const std::string&
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
getCoordinateBases() const
{
  return bases;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
Stokhos::MultiIndex<ordinal_type>
Stokhos::TotalOrderBasis<ordinal_type, value_type, ordering_type>::
getMaxOrders() const
{
  return max_orders;
}
