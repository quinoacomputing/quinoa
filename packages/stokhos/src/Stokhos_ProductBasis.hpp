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

#ifndef STOKHOS_PRODUCTBASIS_HPP
#define STOKHOS_PRODUCTBASIS_HPP

#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {

  /*! 
   * \brief Abstract base class for multivariate orthogonal polynomials
   * generated from tensor products of univariate polynomials.
   */
  /*!
   * * The multivariate polynomials are given by 
   * \f[
   *     \Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)
   * \f]
   * where \f$d\f$ is the dimension of the basis.  This class adds methods 
   * for indexing the multivariate polynomial and getting the coordinate bases.
   */
  template <typename ordinal_type, typename value_type>
  class ProductBasis : 
    public virtual OrthogPolyBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    ProductBasis() {};

    //! Destructor
    virtual ~ProductBasis() {};

    //! Get orders of each coordinate polynomial given an index \c i
    /*!
     * The returned array is of size \f$d\f$, where \f$d\f$ is the dimension of
     * the basis, and entry \f$l\f$ is given by \f$i_l\f$ where
     * \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual Teuchos::Array<ordinal_type> 
    getTerm(ordinal_type i) const = 0;

    //! Get index of the multivariate polynomial given orders of each coordinate
    /*!
     * Given the array \c term storing \f$i_1,\dots,\i_d\f$, returns the index
     * \f$i\f$ such that \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual ordinal_type 
    getIndex(const Teuchos::Array<ordinal_type>& term) const = 0;

    //! Return array of coordinate bases
    /*!
     * Array is of size dimension().
     */
    virtual 
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,
							   value_type> > > 
    getCoordinateBases() const = 0;

  private:

    // Prohibit copying
    ProductBasis(const ProductBasis&);

    // Prohibit Assignment
    ProductBasis& operator=(const ProductBasis& b);

  }; // class ProductBasis

} // Namespace Stokhos

#endif // STOKHOS_PRODUCTBASIS
