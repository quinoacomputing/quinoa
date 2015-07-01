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


module FAztecOO
  use ForTrilinos_enums ,only: FT_AztecOO_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_MultiVector, only: Epetra_MultiVector
  use FEpetra_RowMatrix, only: Epetra_RowMatrix
  use iso_c_binding     ,only: c_int,c_double,c_char
  use foraztecoo
  implicit none
  !private                      ! Hide everything by default
  !public :: AztecOO ! Expose type/constructors/methods

  !> <BR> AztecOO:  An object-oriented wrapper for Aztec.
  !> @{

  !> @brief  AztecOO will solve a linear systems of equations: Ax=b, using Epetra objects and the Aztec solver library, where A is an Epetra_RowMatrix x and b are Epetra_MultiVector objects. 
  type ,extends(universal)                    :: AztecOO !"shell"
  contains
     ! Standard AztecOO solve methods
     !procedure         :: SetAztecOption
     !procedure         :: iterate_current
     !procedure         :: iterate_RowMatrix
     !procedure         :: RecursiveIterate
     !generic :: iterate => iterate_current, iterate_RowMatrix
  end type

contains
  !> <BR> Original C++ prototype:
  !! AztecOO(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_ID_t AztecOO_Create_FromRowMatrix ( CT_Epetra_RowMatrix_ID_t AID,
  ! CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID );
  !> <BR> <BR> ForTrilinos prototype:
  !! AztecOO (Epetra_RowMatrix A, Epetra_MultiVector x, Epetra_MultiVector b);

  !> @name Constructor Function
  !> @{

  type(AztecOO) function AztecOO(A,x,b)
   class(Epetra_RowMatrix) ,intent(in) :: A 
   class(Epetra_MultiVector) ,intent(in) :: x,b 
  end function

  !> @name Constructor Function
  !> @{

  !> @brief Copy Constructor
  type(AztecOO) function AztecOO(this)
    type(AztecOO) ,intent(in) :: this
  end function

  !> @name Standard AztecOO solve methods
  !> @{

  !> @brief AztecOO iteration function.
  !! Iterates on the current problem until MaxIters or Tolerance is reached.
  subroutine iterate(this,MaxIters,tolerance,err) 
    class(AztecOO)   ,intent(in) :: this
    integer(c_int)   ,intent(in) :: MaxIters
    real(c_double)   ,intent(in) :: tolerance
    type(error) ,optional    ,intent(out) :: err
  end subroutine

  !> @name Standard AztecOO solve methods
  !> @{

  !> @brief AztecOO iteration function.
  !! Iterates on the specified matrix and vectors until MaxIters or Toleranceis reached.
  subroutine iterate(this,A,x,b,MaxIters,tolerance,err) 
    class(AztecOO)   ,intent(in) :: this
    class(Epetra_RowMatrix) ,intent(in) :: A
    class(Epetra_MultiVector) ,intent(in) :: x,b
    integer(c_int)   ,intent(in) :: MaxIters
    real(c_double)   ,intent(in) :: tolerance
    type(error) ,optional    ,intent(out) :: err
  end subroutine

  !> @name Special AztecOO solve method
  !> @{

  !> @brief AztecOO iteration functions.
  !! Iterates on the current problem until MaxIters or Tolerance is reached. This one should be suitable for recursive invocations of Aztec.
  subroutine RecursiveIterate(this,MaxIters,tolerance,err) 
    class(AztecOO)   ,intent(in) :: this
    integer(c_int)   ,intent(in) :: MaxIters
    real(c_double)   ,intent(in) :: tolerance
    type(error) ,optional    ,intent(out) :: err
  end subroutine
 
  !> @brief AztecOO option setting function.
  !! Set a specific Aztec option value.
  subroutine SetAztecOption(this,option,value)
   class(AztecOO), intent(in) :: this
   integer(c_int),intent(in) :: option
   integer(c_int),intent(in) :: value
   integer(c_int) ::er
  end subroutine

end module 

