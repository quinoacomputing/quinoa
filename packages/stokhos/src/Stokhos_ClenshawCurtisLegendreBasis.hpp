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

#ifndef STOKHOS_CLENSHAWCURTISLEGENDREBASIS_HPP
#define STOKHOS_CLENSHAWCURTISLEGENDREBASIS_HPP

#include "Stokhos_RecurrenceBasis.hpp"

namespace Stokhos {

  //! Legendre polynomial basis using Clenshaw-Curtis quadrature points
  /*!
   * This is the same as Stokhos::LegendreBasis, but uses Clenshaw-Curtis
   * quadrature points (instead of Gauss-Legendre) for sparse grids only.
   */
  template <typename ordinal_type, typename value_type>
  class ClenshawCurtisLegendreBasis : 
    public RecurrenceBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param normalize whether polynomials should be given unit norm
     */
    ClenshawCurtisLegendreBasis(ordinal_type p, bool normalize = false,
				bool isotropic = false);

    //! Destructor
    ~ClenshawCurtisLegendreBasis();

    //! \name Implementation of Stokhos::OneDOrthogPolyBasis methods
    //@{ 

    /*! 
     * \brief Clone this object with the option of building a higher order
     * basis.
     */
    /*!
     * This method is following the Prototype pattern (see Design Pattern's textbook).
     * The slight variation is that it allows the order of the polynomial to be modified,
     * otherwise an exact copy is formed. The use case for this is creating basis functions
     * for column indices in a spatially varying adaptive refinement context.
     */
    virtual Teuchos::RCP<OneDOrthogPolyBasis<ordinal_type,value_type> > cloneWithOrder(ordinal_type p) const;

    //@}

  protected:

    //! \name Implementation of Stokhos::RecurrenceBasis methods
    //@{ 

    //! Compute recurrence coefficients
    virtual bool
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta,
				  Teuchos::Array<value_type>& gamma) const;

    //@}

    //! Copy constructor with specified order
    ClenshawCurtisLegendreBasis(ordinal_type p, 
				const ClenshawCurtisLegendreBasis& basis);

  private:

    // Prohibit copying
    ClenshawCurtisLegendreBasis(const ClenshawCurtisLegendreBasis&);

    // Prohibit Assignment
    ClenshawCurtisLegendreBasis& operator=(const ClenshawCurtisLegendreBasis& b);

  protected:

    //! Flag determining if expansion is iostropic (same basis in every dim)
    bool isotropic;

  }; // class ClenshawCurtisLegendreBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_ClenshawCurtisLegendreBasisImp.hpp"

#endif
