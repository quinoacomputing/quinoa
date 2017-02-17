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

#ifndef STOKHOS_MP_BLOCK_DIAGONAL_PRECONDITIONER_HPP
#define STOKHOS_MP_BLOCK_DIAGONAL_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"

#include "Stokhos_MPPreconditioner.hpp"
#include "EpetraExt_MultiComm.h"
#include "Epetra_Map.h"
#include "Stokhos_AbstractPreconditionerFactory.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace Stokhos {
    
  /*! 
   * \brief A multi-point preconditioner based on applying the inverse of the
   * diagonal.
   */
  class MPBlockDiagonalPreconditioner : public Stokhos::MPPreconditioner {
      
  public:

    //! Constructor 
    MPBlockDiagonalPreconditioner(
      const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm,
      int num_mp_blocks,
      const Teuchos::RCP<const Epetra_Map>& base_map,
      const Teuchos::RCP<const Epetra_Map>& mp_map,
      const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& prec_factory,
      const Teuchos::RCP<Teuchos::ParameterList>& params);
    
    //! Destructor
    virtual ~MPBlockDiagonalPreconditioner();

    /** \name Stokhos::MPPreconditioner methods */
    //@{

    //! Setup preconditioner
    virtual void 
    setupPreconditioner(
      const Teuchos::RCP<Stokhos::BlockDiagonalOperator>& mp_op, 
      const Epetra_Vector& x);

    //@}

    /** \name Epetra_Operator methods */
    //@{
    
    //! Set to true if the transpose of the operator is requested
    virtual int SetUseTranspose(bool UseTranspose);
    
    /*! 
     * \brief Returns the result of a Epetra_Operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int Apply(const Epetra_MultiVector& Input, 
                      Epetra_MultiVector& Result) const;

    /*! 
     * \brief Returns the result of the inverse of the operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int ApplyInverse(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;
    
    //! Returns an approximate infinity norm of the operator matrix.
    virtual double NormInf() const;
    
    //! Returns a character string describing the operator
    virtual const char* Label () const;
  
    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;
    
    /*! 
     * \brief Returns true if the \e this object can provide an 
     * approximate Inf-norm, false otherwise.
     */
    virtual bool HasNormInf() const;

    /*! 
     * \brief Returns a reference to the Epetra_Comm communicator 
     * associated with this operator.
     */
    virtual const Epetra_Comm & Comm() const;

    /*!
     * \brief Returns the Epetra_Map object associated with the 
     * domain of this matrix operator.
     */
    virtual const Epetra_Map& OperatorDomainMap () const;

    /*! 
     * \brief Returns the Epetra_Map object associated with the 
     * range of this matrix operator.
     */
    virtual const Epetra_Map& OperatorRangeMap () const;

    //@}

  private:
    
    //! Private to prohibit copying
    MPBlockDiagonalPreconditioner(const MPBlockDiagonalPreconditioner&);
    
    //! Private to prohibit copying
    MPBlockDiagonalPreconditioner& operator=(const MPBlockDiagonalPreconditioner&);
    
  protected:
    
    //! Label for operator
    std::string label;

    //! Stores MP parallel communicator
    Teuchos::RCP<const EpetraExt::MultiComm> mp_comm;

    //! Number of mp blocks
    int num_mp_blocks;
    
    //! Stores base map
    Teuchos::RCP<const Epetra_Map> base_map;

    //! Stores MP map
    Teuchos::RCP<const Epetra_Map> mp_map;

    //! Stores factory for building mean preconditioner
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory;

    //! Stores preconditioner for each block
    Teuchos::Array< Teuchos::RCP<Epetra_Operator> > block_precs;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

  }; // class MPBlockDiagonalPreconditioner
  
} // namespace Stokhos

#endif // STOKHOS_MP_BLOCK_DIAGONAL_PRECONDITIONER_HPP
