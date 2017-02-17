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

#include "Ifpack2_OverlapFactor.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_VbrMatrix.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"

//==============================================================================
Ifpack2_OverlapFactor::Ifpack2_OverlapFactor(const Ifpack2_OverlapGraph * OverlapGraph) 
  : Factored_(false),
    Allocated_(false),
    ValuesInitialized_(false),
    OverlapGraph_(OverlapGraph),
    UserMatrix_(0)
{
}
//==============================================================================
Ifpack2_OverlapFactor::Ifpack2_OverlapFactor(const Tpetra_RowMatrix * UserMatrix) 
  : Factored_(false),
    Allocated_(false),
    ValuesInitialized_(false),
    OverlapGraph_(0),
    UserMatrix_(UserMatrix)
{
}
//==============================================================================
Ifpack2_OverlapFactor::Ifpack2_OverlapFactor(const Ifpack2_OverlapFactor & Source) 
  : Factored_(Source.Factored_),
    Allocated_(Source.Allocated_),
    ValuesInitialized_(Source.ValuesInitialized_),
    OverlapGraph_(Source.OverlapGraph_),
    UserMatrix_(Source.UserMatrix_)
{
}
//==============================================================================
int Ifpack2_OverlapFactor::InitValues(const Tpetra_RowMatrix * UserMatrix) {
  

  if (OverlapGraph_!=0) {

    Tpetra_CrsMatrix * CrsMatrix = dynamic_cast<Tpetra_CrsMatrix *>(UserMatrix);
    if (CrsMatrix!=0) 
  if (!Allocated()) EPETRA_CHK_ERR(-1); //Must be allocated
  if (ValuesInitialized()) EPETRA_CHK_ERR(1); // Values already init'ed, warn caller
  
  EPETRA_CHK_ERR(DerivedFactor()); // Call Derived class factorization
  SetValuesInitialized(false);
  SetFactored(true);
  return(0);
}
//==============================================================================
int Ifpack2_OverlapFactor::Factor() {
  
  if (!ValuesInitialized()) EPETRA_CHK_ERR(-1); // Values must be initialized
  if (Factored()) EPETRA_CHK_ERR(1); // Return with a warning that factor already done
  
  EPETRA_CHK_ERR(DerivedFactor()); // Call Derived class factorization
  SetValuesInitialized(false);
  SetFactored(true);
  return(0);
}
