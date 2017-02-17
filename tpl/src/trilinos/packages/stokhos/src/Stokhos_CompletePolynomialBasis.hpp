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

#ifndef STOKHOS_COMPLETEPOLYNOMIALBASIS_HPP
#define STOKHOS_COMPLETEPOLYNOMIALBASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_DerivBasis.hpp"
#include "Stokhos_OneDOrthogPolyBasis.hpp"
#include "Stokhos_ProductBasisUtils.hpp"

namespace Stokhos {

  /*!
   * \brief Multivariate orthogonal polynomial basis generated from a
   * total-order complete-polynomial tensor product of univariate 
   * polynomials.
   */
  /*!
   * The multivariate polynomials are given by 
   * \f[
   *     \Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)
   * \f]
   * where \f$d\f$ is the dimension of the basis and \f$i_1+\dots+ i_d\leq p\f$,
   * where \f$p\f$ is the order of the basis.  The size of the basis is given
   * by \f$(d+p)!/(d!p!)\f$.
   *
   * NOTE:  Currently all coordinate bases must be of the samer order \f$p\f$.
   */
  template <typename ordinal_type, typename value_type>
  class CompletePolynomialBasis : 
    public ProductBasis<ordinal_type,value_type>,
    public DerivBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    /*!
     * \param bases array of 1-D coordinate bases
     * \param sparse_tol tolerance used to drop terms in sparse triple-product
     *                   tensors
     * \param use_old_cijk_alg use old algorithm for computing the sparse
     *                         triple product tensor  (significantly slower,
     *                         but simpler)
     * \param deriv_coeffs direction used to define derivatives for
     *                     derivative product tensors.  Defaults to
     *                     all one's if not supplied.
     */
    CompletePolynomialBasis(
      const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,
 value_type> > >& bases,
      const value_type& sparse_tol = 1.0e-12,
      bool use_old_cijk_alg = false,
      const Teuchos::RCP< Teuchos::Array<value_type> >& deriv_coeffs = Teuchos::null);

    //! Destructor
    virtual ~CompletePolynomialBasis();

    //! \name Implementation of Stokhos::OrthogPolyBasis methods
    //@{

    //! Return order of basis
    ordinal_type order() const;

    //! Return dimension of basis
    ordinal_type dimension() const;

    //! Return total size of basis
    virtual ordinal_type size() const;

    //! Return array storing norm-squared of each basis polynomial
    /*!
     * Entry \f$l\f$ of returned array is given by \f$\langle\Psi_l^2\rangle\f$
     * for \f$l=0,\dots,P\f$ where \f$P\f$ is size()-1.
     */
    virtual const Teuchos::Array<value_type>& norm_squared() const;

    //! Return norm squared of basis polynomial \c i.
    virtual const value_type& norm_squared(ordinal_type i) const;

    //! Compute triple product tensor
    /*!
     * The \f$(i,j,k)\f$ entry of the tensor \f$C_{ijk}\f$ is given by
     * \f$C_{ijk} = \langle\Psi_i\Psi_j\Psi_k\rangle\f$ where \f$\Psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j,k=0,\dots,P\f$ where
     * \f$P\f$ is size()-1.
     */
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensor() const;

    //! Compute linear triple product tensor where k = 0,1,..,d
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeLinearTripleProductTensor() const;

    //! Evaluate basis polynomial \c i at zero
    virtual value_type evaluateZero(ordinal_type i) const;

    //! Evaluate basis polynomials at given point \c point
    /*!
     * Size of returned array is given by size(), and coefficients are
     * ordered from order 0 up to size size()-1.
     */
    virtual void evaluateBases(
      const Teuchos::ArrayView<const value_type>& point,
      Teuchos::Array<value_type>& basis_vals) const;

    //! Print basis to stream \c os
    virtual void print(std::ostream& os) const;

    //! Return string name of basis
    virtual const std::string& getName() const;

    //@}

    //! \name Implementation of Stokhos::ProductBasis methods
    //@{

    //! Get orders of each coordinate polynomial given an index \c i
    /*!
     * The returned array is of size \f$d\f$, where \f$d\f$ is the dimension of
     * the basis, and entry \f$l\f$ is given by \f$i_l\f$ where
     * \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual const MultiIndex<ordinal_type>& term(ordinal_type i) const;

    //! Get index of the multivariate polynomial given orders of each coordinate
    /*!
     * Given the array \c term storing \f$i_1,\dots,\i_d\f$, returns the index
     * \f$i\f$ such that \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual ordinal_type index(const MultiIndex<ordinal_type>& term) const;

    //! Return coordinate bases
    /*!
     * Array is of size dimension().
     */
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, 
							   value_type> > > 
    getCoordinateBases() const;

    //! Return maximum order allowable for each coordinate basis
    virtual MultiIndex<ordinal_type> getMaxOrders() const;

    //@}

    //! \name Implementation of Stokhos::DerivBasis methods
    //@{

    /*! 
     * \brief Compute triple product tensor 
     * \f$D_{ijk} = \langle\Psi_i\Psi_j D_v\Psi_k\rangle\f$ where 
     * \f$D_v\Psi_k\f$ represents the derivative of \f$\Psi_k\f$ in the 
     * direction \f$v\f$.
     */
    /*!
     * The definition of \f$v\f$ is defined by the \c deriv_coeffs 
     * constructor argument.
     */
    virtual 
    Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > 
    computeDerivTripleProductTensor(
      const Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >& Bij,
      const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk
      ) const;

    /*! 
     * \brief Compute double product tensor 
     * \f$B_{ij} = \langle \Psi_i D_v\Psi_j\rangle\f$ where \f$D_v\Psi_j\f$
     * represents the derivative of \f$\Psi_j\f$ in the direction \f$v\f$.
     */
    /*!
     * The definition of \f$v\f$ is defined by  the \c deriv_coeffs 
     * constructor argument.
     */
    virtual 
    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > 
    computeDerivDoubleProductTensor() const;

    //@}

  protected:

    //! Compute triple product tensor using old algorithm
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensorOld(ordinal_type order) const;

    //! Compute triple product tensor using new algorithm
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensorNew(ordinal_type order) const;

  private:

    // Prohibit copying
    CompletePolynomialBasis(const CompletePolynomialBasis&);

    // Prohibit Assignment
    CompletePolynomialBasis& operator=(const CompletePolynomialBasis& b);
    
  protected:

    typedef Stokhos::CompletePolynomialBasisUtils<ordinal_type,value_type> CPBUtils;

    //! Name of basis
    std::string name;

    //! Total order of basis
    ordinal_type p;

    //! Total dimension of basis
    ordinal_type d;

    //! Total size of basis
    ordinal_type sz;

    //! Array of bases
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > > bases;

    //! Array storing order of each basis
    Teuchos::Array<ordinal_type> basis_orders;

    //! Tolerance for computing sparse Cijk
    value_type sparse_tol;

    //! Use old algorithm for computing Cijk
    bool use_old_cijk_alg;

    //! Coefficients for derivative
    Teuchos::RCP< Teuchos::Array<value_type> > deriv_coeffs;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! 2-D array of basis terms
    Teuchos::Array< MultiIndex<ordinal_type> > terms;

    //! Number of terms up to each order
    Teuchos::Array<ordinal_type> num_terms;

    //! Temporary array used in basis evaluation
    mutable Teuchos::Array< Teuchos::Array<value_type> > basis_eval_tmp;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;

  }; // class CompletePolynomialBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_CompletePolynomialBasisImp.hpp"

#endif
