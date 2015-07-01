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

module FEpetra_Vector
  use ForTrilinos_enums   ,only: FT_Epetra_MultiVector_ID_t,FT_Epetra_Vector_ID_t,&
                                 FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t,FT_boolean_t 
  use ForTrilinos_table_man
  use ForTrilinos_error
  use ForTrilinos_universal
  use FEpetra_MultiVector ,only: Epetra_MultiVector
  use FEpetra_BlockMap    !,only: Epetra_BlockMap !use to circumvent reported compiler bug
  use iso_c_binding       ,only: c_int
  use forepetra
  implicit none
  private                                    ! Hide everything by default
  public :: Epetra_Vector !,Epetra_MultiVector ! Expose type/constructors/methods

  type ,extends(Epetra_MultiVector)      :: Epetra_Vector !"shell"
  contains
     !Developers only
     procedure         :: invalidate_id => invalidate_EpetraVector_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraVector
     procedure         :: get_EpetraVector_ID 
     procedure ,nopass :: alias_EpetraVector_ID
     procedure         :: generalize 
     !Post-construction modfication routines
     procedure ,private:: ReplaceGlobalValues_NoOffset
     procedure ,private:: ReplaceGlobalValues_BlockPos
     generic :: ReplaceGlobalValues => ReplaceGlobalValues_NoOffset,ReplaceGlobalValues_BlockPos
     procedure ,private:: ReplaceMyValues_NoOffset
     procedure ,private:: ReplaceMyValues_BlockPos
     generic :: ReplaceMyValues => ReplaceMyValues_NoOffset,ReplaceMyValues_BlockPos
     procedure         :: SumIntoGlobalValues_NoOffset
     procedure         :: SumIntoGlobalValues_BlockPos
     generic :: SumIntoGlobalValues => SumIntoGlobalValues_NoOffset,SumIntoGlobalValues_BlockPos
     procedure         :: SumIntoMyValues_NoOffset
     procedure         :: SumIntoMyValues_BlockPos
     generic :: SumIntoMyValues => SumIntoMyValues_NoOffset,SumIntoMyValues_BlockPos
     ! Extraction methods
     procedure         :: ExtractCopy_EpetraVector
     generic :: ExtractCopy => ExtractCopy_EpetraVector
     !overloaded operators
     procedure         :: get_element_EpetraVector
     generic :: get_Element => get_element_EpetraVector 
  end type

 
contains


  ! Original C++ prototype:
  ! Epetra_Vector(const Epetra_BlockMap& Map, bool zeroOut = true);
  ! CTrilinos prototype:
  ! CT_Epetra_Vector_ID_t Epetra_Vector_Create ( CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut );

  !> @name Constructor Functions
  !! @{
 
  !> <BR> Epetra_Vector constructor conformal to a BlockMap, optionally zero the newly created vector.
  type(Epetra_Vector) function Epetra_Vector(BlockMap,zero_initial)
    type(Epetra_BlockMap) ,intent(in) :: BlockMap
    logical ,optional      ,intent(in) :: zero_initial
  end function
  
  !> @name Constructor Functions
  !! @{
 
  !> <BR> Epetra_Vector constructor from a user supplied vector
  type(Epetra_Vector) function Epetra_Vector(CV,BlockMap,V)
    type(Epetra_BlockMap) ,intent(in) :: BlockMap
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    real(c_double),dimension(:) :: V
  end function

  ! Original C++ prototype:
  ! Epetra_Vector(const Epetra_Vector& Source);
  ! CTrilinos prototype:
  ! CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( CT_Epetra_Vector_ID_t SourceID );

  !> @name Constructor Functions
  !! @{
 
  !> <BR> Epetra_Vector copy constructor
  type(Epetra_Vector) function Epetra_Vector(this)
    type(Epetra_Vector) ,intent(in) :: this 
  end function


  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces entries with a list of indexed values, indices in global space
  subroutine ReplaceGlobalValues(this,values,indices,err)
    type(Epetra_Vector), intent(in) :: this
    real(c_double),dimension(:),intent(in) :: values
    integer(c_int),dimension(:),intent(in) :: indices 
    type(error),optional,intent(out) :: err
  end subroutine

  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces entries with a list of indexed values, indices in global space
  subroutine ReplaceGlobalValues(this,BlockOffset,values,indices,err)
    type(Epetra_Vector), intent(in) :: this
    integer(c_int)       ,intent(in) :: BlockOffset
    real(c_double),dimension(:),intent(in) :: values
    integer(c_int),dimension(:),intent(in) :: indices 
    type(error),optional,intent(out) :: err
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces entries with a list of indexed values, indices in local space
  subroutine ReplaceMyValues(this,values,indices,err)
    type(Epetra_Vector), intent(in) :: this
    real(c_double),dimension(:),intent(in) :: values
    integer(c_int),dimension(:),intent(in) :: indices 
    integer(c_int),dimension(size(indices)):: indices_c
    type(error),optional,intent(out) :: err
  end subroutine

  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces entries with a list of indexed values, indices in local space
  subroutine ReplaceMyValues(this,BlockOffset,values,indices,err)
    type(Epetra_Vector), intent(in) :: this
    integer(c_int)       ,intent(in) :: BlockOffset
    real(c_double),dimension(:),intent(in) :: values
    integer(c_int),dimension(:),intent(in) :: indices 
    integer(c_int),dimension(size(indices)):: indices_c
    type(error),optional,intent(out) :: err
  end subroutine

  !> @name Post-construction modification routines
  !! @{

  !> <BR> Adds entries with list of indexed values, indices in global space
  subroutine SumIntoGlobalValues(this,values,indices,err)
    type(Epetra_Vector), intent(in) :: this
    real(c_double),dimension(:),intent(in) :: values
    integer(c_int),dimension(:),intent(in) :: indices 
    type(error),optional,intent(out) :: err
    integer(c_int)                      :: error_out
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR> Adds entries with list of indexed values, indices in global space
  subroutine SumIntoGlobalValues(this,BlockOffset,values,indices,err)
    type(Epetra_Vector), intent(in) :: this
    integer(c_int)       ,intent(in) :: BlockOffset
    real(c_double),dimension(:),intent(in) :: values
    integer(c_int),dimension(:),intent(in) :: indices 
    type(error),optional,intent(out) :: err
  end subroutine

  !> @name Post-construction modification routines
  !! @{

  !> <BR> Adds entries with list of indexed values, indices in local space
  subroutine SumIntoMyValues(this,values,indices,err)
    type(Epetra_Vector), intent(in) :: this
    real(c_double),dimension(:),intent(in) :: values
    integer(c_int),dimension(:),intent(in) :: indices 
    integer(c_int),dimension(size(indices)):: indices_c   
    type(error),optional,intent(out) :: err
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR> Adds entries with list of indexed values, indices in local space
  subroutine SumIntoMyValues(this,BlockOffset,values,indices,err)
    type(Epetra_Vector), intent(in) :: this
    integer(c_int)       ,intent(in) :: BlockOffset
    real(c_double),dimension(:),intent(in) :: values
    integer(c_int),dimension(:),intent(in) :: indices 
    type(error),optional,intent(out) :: err
  end subroutine

  !> @name Extraction methods
  !! @{

  !> <BR> Copies vector contents into target
  function ExtractCopy(this,err) result(ExtractCopy_out)
    type(Epetra_Vector), intent(in) :: this
    real(c_double), dimension(:), allocatable::  ExtractCopy_out 
    type(error),optional,intent(out) :: err
  end function 
 
  !> @name Extraction methods
  !! @{

  !> <BR> Extract value from an entry in  vector
  real(c_double) function get_element_EpetraVector(this,index)
    type(Epetra_Vector), intent(in) :: this
    integer(c_int), intent(in) :: index
  end function
   
end module 

