/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_REORDERING_HPP
#define IFPACK2_REORDERING_HPP

#include "Ifpack2_ConfigDefs.hpp"

namespace Teuchos {
  class ParameterList;
}
class Tpetra_MultiVector;
class Ifpack2_Graph;
class Tpetra_RowMatrix;

//! Ifpack2_Reordering: basic class for reordering for a Ifpack2_Graph object.
/*!
Class Ifpack2_Reordering is a pure virtual class that defines the 
structure of all Ifpack2 reordering.

The Ifpack2_Graph object is used \e only by method Compute().

A typical code reads as follows (using for instance RCM reordering):
\code
#include "Ifpack2_Reordering.hpp"
#include "Ifpack2_RCMReordering.hpp"
#include "Ifpack2_Graph.hpp"
// A is an Tpetra_RowMatrix pointer.
// Need to create a graph (which is a simple wrapper)
// This required include file Ifpack2_Graph_Tpetra_RowMatrix.hpp
Ifpack2_Graph_Tpetra_RowMatrix Graph(A);

// Construct the object
Ifpack2_RCMReordering Reorder(Graph);
// Compute the reordering.
IFPACK2_CHK_ERR(Reorder.Compute());
// Prints out some information
cout << Reorder;
\endcode

<P>An Ifpack2_Reordering object is a tool used by class Ifpack2_Preconditioner,
to reorder the localized matrix (with or without overlap). As its
basic usage is for localized matrices, this class takes care of
reordering the \e local rows only. It is also supposed that the
input graph contains no singletons. This is not a limitation, as
class Ifpack2_AdditiveSchwarz will filter the graph using
Ifpack2_SingletonFilter before using reordering.

<P>If TIFPACK is configure with Teuchos support, method SetParameters()
should be adopted. Otherwise, users can set parameters (one at-a-time),
using methods SetParameter(), for integers and doubles.

<P>Ifpack2_Preconditioner objects overload the << operator. Derived
classes should specify a Print() method, that will be used in
operator <<.

\author Michael Heroux, SNL 9214.

\date Last modified: Oct-04.
*/

class Ifpack2_Reordering {

public:

  //! Destructor.
  virtual ~Ifpack2_Reordering() {};
  
  //! Sets integer parameters `Name'.
  virtual int SetParameter(const string Name, const int Value) = 0;
 
  //! Sets double parameters `Name'.
  virtual int SetParameter(const string Name, const double Value) = 0;

  //! Sets all parameters.
  virtual int SetParameters(Teuchos::ParameterList& List) = 0;

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Ifpack2_Graph& Graph) = 0;

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Tpetra_RowMatrix& Matrix) = 0;

  //! Returns \c true is the reordering object has been successfully initialized, false otherwise.
  virtual bool IsComputed() const = 0;

  //! Returns the reordered index of row \c i.
  virtual int Reorder(const int i) const = 0;

  //! Returns the inverse reordered index of row \c i.
  virtual int InvReorder(const int i) const = 0;

  //! Applies reordering to multivector Xorig, whose local length equals the number of local rows, stores reordered vector in X.
  virtual int P(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xorig,
		Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X) const = 0;

  //! Applies inverse reordering to multivector Xorig, whose local length equals the number of local rows, stores inverse reordered vector in X.
  virtual int Pinv(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xorig,
		   Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X) const = 0;

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual ostream& Print(std::ostream& os) const = 0;

}; 

inline ostream& operator<<(ostream& os, const Ifpack2_Reordering& obj)
{
  return(obj.Print(os));
}

#endif
