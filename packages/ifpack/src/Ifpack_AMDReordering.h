/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_AMDREORDERING_H
#define IFPACK_AMDREORDERING_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Reordering.h"

namespace Teuchos {
  class ParameterList;
}
class Ifpack_Graph;
class Epetra_MultiVector;
class Epetra_RowMatrix;

//! Ifpack_AMDReordering: approximate minimum degree reordering.

class Ifpack_AMDReordering : public Ifpack_Reordering {

public:

  //! Constructor for Ifpack_Graph's.
  Ifpack_AMDReordering();

  //! Copy Constructor.
  Ifpack_AMDReordering(const Ifpack_AMDReordering& RHS);

  //! Assignment operator.
  Ifpack_AMDReordering& operator=(const Ifpack_AMDReordering& RHS);

  //! Destructor.
  virtual ~Ifpack_AMDReordering() {};
  
  //! Sets integer parameters `Name'.
  int SetParameter(const string Name, const int Value);

  //! Sets double parameters `Name'.
  int SetParameter(const string Name, const double Value);
  
  //! Sets all parameters.
  int SetParameters(Teuchos::ParameterList& List);

  //! Computes all it is necessary to initialize the reordering object.
  int Compute(const Ifpack_Graph& Graph);

  //! Computes all it is necessary to initialize the reordering object.
  int Compute(const Epetra_RowMatrix& Matrix);

  //! Returns \c true is the reordering object has been successfully initialized, false otherwise.
  bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Returns the reordered index of row \c i.
  inline int Reorder(const int i) const;

  //! Returns the inverse reordered index of row \c i.
  inline int InvReorder(const int i) const;

  //! Applies reordering to multivector X, whose local length equals the number of local rows.
  int P(const Epetra_MultiVector& Xorig,
	Epetra_MultiVector& Xreord) const;

  //! Applies inverse reordering to multivector X, whose local length equals the number of local rows.
  int Pinv(const Epetra_MultiVector& Xorig,
	   Epetra_MultiVector& Xinvreord) const;

  
  //! Prints basic information on iostream. This function is used by operator<<.
  ostream& Print(std::ostream& os) const;

  //! Returns the number of local rows.
  int NumMyRows() const 
  {
    return(NumMyRows_);
  }

private:
  //! Number of local rows in the graph.
  int NumMyRows_;
  //! If \c true, the reordering has been successfully computed.
  bool IsComputed_;
  //! Contains the reordering.
  std::vector<int> Reorder_;
  //! Contains the inverse reordering.
  std::vector<int> InvReorder_;
}; 

#endif
