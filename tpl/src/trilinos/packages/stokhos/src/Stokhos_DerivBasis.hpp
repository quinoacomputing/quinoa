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

#ifndef STOKHOS_DERIVBASIS_HPP
#define STOKHOS_DERIVBASIS_HPP

#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {
  
  /*! 
   * \brief Abstract base class for multivariate orthogonal polynomials
   * that support computing double and triple products involving derivatives
   * of the basis polynomials.
   */
  template <typename ordinal_type, typename value_type>
  class DerivBasis : 
    public virtual OrthogPolyBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    DerivBasis() {};

    //! Destructor
    virtual ~DerivBasis() {};

    /*! 
     * \brief Compute triple product tensor 
     * \f$D_{ijk} = \langle\Psi_i\Psi_j D_v\Psi_k\rangle\f$ where 
     * \f$D_v\Psi_k\f$ represents the derivative of \f$\Psi_k\f$ in the 
     * direction \f$v\f$.
     */
    /*!
     * The definition of \f$v\f$ is defined by the derived class implementation.
     */
    virtual 
    Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > 
    computeDerivTripleProductTensor(
      const Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >& Bij,
      const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk) const = 0;

    /*! 
     * \brief Compute double product tensor 
     * \f$B_{ij} = \langle \Psi_i D_v\Psi_j\rangle\f$ where \f$D_v\Psi_j\f$
     * represents the derivative of \f$\Psi_j\f$ in the direction \f$v\f$.
     */
    /*!
     * The definition of \f$v\f$ is defined by the derived class implementation.
     */
    virtual 
    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > 
    computeDerivDoubleProductTensor() const = 0;

  private:

    // Prohibit copying
    DerivBasis(const DerivBasis&);

    // Prohibit Assignment
    DerivBasis& operator=(const DerivBasis& b);

  }; // class DerivBasis

} // Namespace Stokhos

#endif // STOKHOS_DERIVBASIS
