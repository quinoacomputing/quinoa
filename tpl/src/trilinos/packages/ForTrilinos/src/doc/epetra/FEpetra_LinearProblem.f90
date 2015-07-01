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


module FEpetra_LinearProblem
  use ForTrilinos_enums !,only: FT_Epetra_LinearProblem_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_universal, only: universal
  use ForTrilinos_table_man
  use ForTrilinos_error
  use FEpetra_MultiVector, only: Epetra_MultiVector
  use FEpetra_RowMatrix, only: Epetra_RowMatrix
  !use FEpetra_Operator, only: Epetra_Operator
  use iso_c_binding      ,only: c_int
  use forepetra
  implicit none
  private                     ! Hide everything by default
  public :: Epetra_LinearProblem ! Expose type/constructors/methods

  !> <BR> Epetra_LinearProblem:  The Epetra Linear Problem Class. 
  !> @{

  !> @brief The Epetra_LinearProblem class is a wrapper that encapsulates the general information needed for solving a linear system of equations.  !! Currently it accepts a Epetra_RowMatrix, initial guess and RHS and returns the solution. 
  type :: Epetra_LinearProblem    !, extends(universal)      :: Epetra_LinearProblem 
  contains
  end type

 
contains


  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_LinearProblem default construtor.
  !> @brief Creates an empty Epetra_LinearProblem object. The operator A, left hand side X and right hand side B  must be set using the SetOperator(), SetLhs() and SetRhs() methods. 
  type(Epetra_LinearProblem) function Epetra_LinearProblem()
  end function

  
  
  !> @name Constructor Functions
  !! @{

  !> <BR>  Epetra_LinearProblem  Constructor passing the operator as a matrix.
  !> @brief Creates an Epetra_LinearProblem where the operator A is passed as a matrix.
  type(Epetra_LinearProblem) function Epetra_LinarProblem(A,X,B)
    class(Epetra_RowMatrix), intent(in) :: A  &
         !< In The operator matrix

    class(Epetra_MultiVector), intent(in) :: X &
    !< In The left hand side X

    class(Epetra_MultiVector), intent(in) :: B &
    !< In The right hand side B.

  end function

  
  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_LinearProblem Copy constructor.
  type(Epetra_LinearProblem) function Epetra_LinearProblem(this)
    type(Epetra_LinearProblem) ,intent(in) :: this
  end function


  !> @name Integrity Check Methods
  !! @{

  !> <BR>   Sanity tests on input problem.
  integer(c_int) function CheckInput(this)
    class(Epetra_LinearProblem), intent(in) :: this
  end function
 
  !> @name Set Methods
  !! @{

  !> <BR>   Set assertion of symmetry for current problem .
  subroutine AssertSymmetric(this)
    class(Epetra_LinearProblem), intent(in) :: this
  end subroutine

  !> @name Set Methods
  !! @{

  !> <BR>   Set problem difficulty level.
  subroutine SetPDL(this,PDL)
    class(Epetra_LinearProblem), intent(in) :: this
    integer(FT_ProblemDifficultyLevel_E_t), intent(in) :: PDL
  end subroutine SetPDL

  !> @name Set Methods
  !! @{

  !> <BR>   Set operator matrix A.
 subroutine SetOperator_Matrix(this,A)
   class(Epetra_LinearProblem), intent(in) :: this
   class(Epetra_RowMatrix), intent(in) :: A
 end subroutine

  !> @name Set Methods
  !! @{

  !> <BR>   Set left hand side X. 
 subroutine SetLHS(this,X)
   class(Epetra_LinearProblem), intent(in) :: this
   class(Epetra_MultiVector), intent(in) :: X
 end subroutine
  
  !> @name Set Methods
  !! @{

  !> <BR>  Set right  hand side B.
 subroutine SetRHS(this,B)
   class(Epetra_LinearProblem), intent(in) :: this
   class(Epetra_MultiVector), intent(in) :: B
 end subroutine


end module 

