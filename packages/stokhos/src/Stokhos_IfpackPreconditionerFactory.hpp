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

#ifndef STOKHOS_IFPACK_PRECONDITIONER_FACTORY_HPP
#define STOKHOS_IFPACK_PRECONDITIONER_FACTORY_HPP

#include "Teuchos_ParameterList.hpp"

#include "Stokhos_AbstractPreconditionerFactory.hpp"

namespace Stokhos {

  //! A factory for building Ifpack preconditioners
  class IfpackPreconditionerFactory : 
    public Stokhos::AbstractPreconditionerFactory {
  public:

    //! Constructor
    IfpackPreconditionerFactory(const Teuchos::RCP<Teuchos::ParameterList>& p);

    //! Destructor
    virtual ~IfpackPreconditionerFactory() {}

    //! Compute preconditioner
    virtual Teuchos::RCP<Epetra_Operator> 
    compute(const Teuchos::RCP<Epetra_Operator>& op,
	    bool compute_prec = true);

    //! Recompute preconditioner operator for a new matrix
    virtual void
    recompute(const Teuchos::RCP<Epetra_Operator>& op,
	      const Teuchos::RCP<Epetra_Operator>& prec);

  protected:

    //! Preconditioner parameters
    Teuchos::RCP<Teuchos::ParameterList> precParams;

  }; // class IfpackPreconditionerFactory

} // namespace Stokhos

#endif // STOKHOS_IFPACK_PRECONDITIONER_FACTORY_HPP
