!*********************************************************************
! ForTrilinos: Object-Oriented Fortran 2003 interface to Trilinos
!                Copyright 2010 Sandia Corporation
!
! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
! the U.S. Government retains certain rights in this software.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the Corporation nor the names of the
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Questions? Contact Karla Morris  (knmorri@sandia.gov) or
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************


module FEpetra_RowMatrix
  use ForTrilinos_universal,only : universal
  use ForTrilinos_enums !,only: FT_Epetra_RowMatrix_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_error
  use FEpetra_Map, only: Epetra_Map
  use FEpetra_MultiVector, only: Epetra_MultiVector
  use FEpetra_Vector, only: Epetra_Vector
  use forepetra
  implicit none
  !private               ! Hide everything by default
  !public :: Epetra_RowMatrix ! Expose type/methods

  !> <BR> Epetra_RowMatrix: An abstract class for using real-valued double-precision row matrices.
  !! @{

  !> @brief The Epetra_RowMatrix class is an abstract class (specifies interface only) that enable the use of real-valued double-precision sparse matrices where matrix entries are intended for row access.  It is currently implemented by the Epetra_CrsMatrix class.

  type ,abstract :: Epetra_RowMatrix !,extends(universal) :: Epetra_RowMatrix
  contains
    !Matrix data extraction routines
    !procedure(NumMyRowEntries_interface),deferred :: NumMyRowEntries
    !procedure(MaxNumEntries_interface)  ,deferred :: MaxNumEntries
    ! Computational Methods
    !procedure(Multiply_interface) ,deferred :: Multiply
    !Attribute access functions
    !procedure(RowMatrixRowMap_interface),deferred :: RowMatrixRowMap
  end type

  abstract interface
    !> @name Matrix Data Extraction Routines
    !! @{
 
    !> @brief Returns the number of nonzero entries in MyRow.
    integer(c_int) function NumMyRowEntries(this,MyRow)
      use iso_c_binding, only : c_int
      import:: Epetra_RowMatrix
      class(Epetra_RowMatrix), intent(in) :: this 
      integer(c_int),          intent(in) :: MyRow &
      !<  Local row
    end function
    !> @name Matrix Data Extraction Routines
    !! @{
 
    !> @brief Returns the maximum of NumMyRowEntries() over all rows.
    integer(c_int) function MaxNumEntries(this)
      use iso_c_binding, only : c_int
      import:: Epetra_RowMatrix
      class(Epetra_RowMatrix), intent(in) :: this
    end function
    !> @name Mathematical Functions
    !! @{
 
    !> @brief Returns the result of a Epetra_RowMatrix multiplied by a Epetra_MultiVector X in Y.
     subroutine Multiply(this,TransA,x,y,err)
      use iso_c_binding, only: c_int
      import :: Epetra_RowMatrix,Epetra_MultiVector,error
      class(Epetra_RowMatrix), intent(in) :: this
      logical, intent(in) :: TransA &
      !< If true, multiply by the transpose of matrix, otherwise just use matrix.
      class(Epetra_MultiVector), intent(in) :: x &
      !< A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
      class(Epetra_MultiVector), intent(inout) :: y &
      !< A Epetra_MultiVector of dimension NumVectorscontaining result.
      type(error), optional,intent(inout) :: err &
      !< Returns  error information.
     end subroutine
    !> @name Attribute access Function
    !! @{
 
    !> @brief Returns the Epetra_Map object associated with the rows of this matrix.
    function RowMatrixRowMap(this) 
     import:: Epetra_RowMatrix,Epetra_Map
     class(Epetra_RowMatrix), intent(in) :: this
     type(Epetra_Map) :: RowMatrixRowMap
    end function 
  end interface

end module 
