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

#ifndef STOKHOS_EPETRAVECTORORTHOGPOLY_HPP
#define STOKHOS_EPETRAVECTORORTHOGPOLY_HPP

#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_ProductEpetraVector.hpp"

namespace Stokhos {

  /*! 
   * \brief A container class storing an orthogonal polynomial whose
   * coefficients are vectors, operators, or in general any type that 
   * would have an expensive copy constructor.  
   */
  class EpetraVectorOrthogPoly : 
    public VectorOrthogPoly<Epetra_Vector>,
    public ProductEpetraVector {
  public:

    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    //! Constructor with no basis
    /*!
     * Use with care!  Generally you will want to call reset() before using
     * any of the methods on this class.
     */
    EpetraVectorOrthogPoly();

    /*! 
     * \brief Create a polynomial for basis \c basis with empty 
     * coefficients
     */
    EpetraVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * generated from the supplied map.
     */
    EpetraVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * given by a created block vector
     */
    EpetraVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * given by the supplied block vector.
     */
    EpetraVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm,
      Epetra_DataAccess CV,
      const Epetra_Vector& block_vector);
    
    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraVectorOrthogPoly(const EpetraVectorOrthogPoly& v);

    //! Destructor
    virtual ~EpetraVectorOrthogPoly();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraVectorOrthogPoly& operator=(const EpetraVectorOrthogPoly& v);

    //! Reset to a new basis
    /*!
     * This resizes array to fit new basis.
     */
    void reset(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    //! Reset to a new basis
    /*!
     * This resizes array to fit new basis.
     */
    void reset(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    //! Compute mean
    void computeMean(Epetra_Vector& v) const;

    //! Compute variance
    void computeVariance(Epetra_Vector& v) const;

    //! Compute standard deviation
    void computeStandardDeviation(Epetra_Vector& v) const;

  }; // class EpetraVectorOrthogPoly

} // end namespace Stokhos

#endif  // STOKHOS_EPETRAVECTORORTHOGPOLY_HPP
