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

#ifndef STOKHOS_COMPLETEPOLYNOMIALBASIS_HPP
#define STOKHOS_COMPLETEPOLYNOMIALBASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_DerivBasis.hpp"
#include "Stokhos_OneDOrthogPolyBasis.hpp"

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
      const value_type& sparse_tol = 1.0e-15,
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
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is size()-1 and \f$k=0,\dots,p\f$ where \f$p\f$
     * is the supplied \c order.
     */
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensor(ordinal_type order) const;

    //! Evaluate basis polynomial \c i at zero
    virtual value_type evaluateZero(ordinal_type i) const;

    //! Evaluate basis polynomials at given point \c point
    /*!
     * Size of returned array is given by size(), and coefficients are
     * ordered from order 0 up to size size()-1.
     */
    virtual void evaluateBases(const Teuchos::Array<value_type>& point,
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
    virtual Teuchos::Array<ordinal_type> getTerm(ordinal_type i) const;

    //! Get index of the multivariate polynomial given orders of each coordinate
    /*!
     * Given the array \c term storing \f$i_1,\dots,\i_d\f$, returns the index
     * \f$i\f$ such that \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual ordinal_type 
    getIndex(const Teuchos::Array<ordinal_type>& term) const;

    //! Return coordinate bases
    /*!
     * Array is of size dimension().
     */
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, 
							   value_type> > > 
    getCoordinateBases() const;

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

    /*!
     * \brief Compute the 2-D array of basis terms which maps a basis index
     * into the orders for each basis dimension
     */
    void compute_terms();

    /*!
     * \brief Compute basis index given the orders for each basis
     * dimension.
     */
    ordinal_type compute_index(const Teuchos::Array<ordinal_type>& terms) const;

  private:

    // Prohibit copying
    CompletePolynomialBasis(const CompletePolynomialBasis&);

    // Prohibit Assignment
    CompletePolynomialBasis& operator=(const CompletePolynomialBasis& b);
    
  protected:

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
    Teuchos::Array< Teuchos::Array<ordinal_type> > terms;

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
