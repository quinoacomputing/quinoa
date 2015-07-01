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
  private               ! Hide everything by default
  public :: Epetra_RowMatrix ! Expose type/methods

  type ,abstract ,extends(universal) :: Epetra_RowMatrix
   private
   type(FT_Epetra_RowMatrix_ID_t)         :: RowMatrix_id 
  contains
    !Matrix data extraction routines
    procedure(NumMyRowEntries_interface),deferred :: NumMyRowEntries
    procedure(MaxNumEntries_interface)  ,deferred :: MaxNumEntries
    ! Computational Methods
    procedure(Multiply_interface) ,deferred :: Multiply
    !Attribute access functions
    procedure(RowMatrixRowMap_interface),deferred :: RowMatrixRowMap
    ! ForTriliinos developers only -- not for invocation from end-user application code:
    procedure                          :: invalidate_EpetraRowMatrix_ID
    procedure                          :: ctrilinos_delete_EpetraRowMatrix
    procedure                          :: get_EpetraRowMatrix_ID
    procedure                          :: set_EpetraRowMatrix_ID
    procedure                 ,nopass  :: alias_EpetraRowMatrix_ID
    procedure ,non_overridable         :: generalize_EpetraRowMatrix
  end type
  
  abstract interface
    integer(c_int) function NumMyRowEntries_interface(this,MyRow)
      use iso_c_binding, only : c_int
      import:: Epetra_RowMatrix
      class(Epetra_RowMatrix), intent(in) :: this
      integer(c_int),          intent(in) :: MyRow
    end function
    integer(c_int) function MaxNumEntries_interface(this)
      use iso_c_binding, only : c_int
      import:: Epetra_RowMatrix
      class(Epetra_RowMatrix), intent(in) :: this
    end function
     subroutine Multiply_interface(this,TransA,x,y,err)
      use iso_c_binding, only: c_int
      import :: Epetra_RowMatrix,Epetra_MultiVector,error
      class(Epetra_RowMatrix), intent(in) :: this
      logical, intent(in) :: TransA
      class(Epetra_MultiVector), intent(in) :: x
      class(Epetra_MultiVector), intent(inout) :: y
      type(error), optional,intent(inout) :: err
     end subroutine
    function RowMatrixRowMap_interface(this) 
     import:: Epetra_RowMatrix,Epetra_Map
     class(Epetra_RowMatrix), intent(in) :: this
     type(Epetra_Map) :: RowMatrixRowMap_interface
    end function 
  end interface

  contains
  
  type(FT_Epetra_RowMatrix_ID_t) function get_EpetraRowMatrix_ID(this)
    class(Epetra_RowMatrix) ,intent(in) :: this
    get_EpetraRowMatrix_ID = this%RowMatrix_id
  end function
  
  subroutine set_EpetraRowMatrix_ID(this,id)
    class(Epetra_RowMatrix)        ,intent(inout) :: this
    type(FT_Epetra_RowMatrix_ID_t) ,intent(in)    :: id 
    this%RowMatrix_id=id
  end subroutine 
  
  type(FT_Epetra_RowMatrix_ID_t) function alias_EpetraRowMatrix_ID(generic_id)
    use iso_c_binding, only : c_loc,c_int
    use ForTrilinos_table_man
    use ForTrilinos_enums
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_RowMatrix_ID),stat=status)
    ierr=error(status,'FEpetra_RowMatrix:alias_EpetraRowMatrix_ID')
    call ierr%check_success()
    alias_EpetraRowMatrix_ID=degeneralize_EpetraRowMatrix(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize_EpetraRowMatrix(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_RowMatrix) ,intent(in) ,target :: this
   generalize_EpetraRowMatrix = generalize_all( c_loc(this%RowMatrix_id) )
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_RowMatrix) ,intent(in) ,target :: this
   ! generalize_EpetraRowMatrix = Epetra_RowMatrix_Generalize ( this%RowMatrix_id )
   ! ____ Use for CTrilinos function implementation ______
  end function
  
  type(FT_Epetra_RowMatrix_ID_t) function degeneralize_EpetraRowMatrix(generic_id) bind(C)
    !use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_RowMatrix_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)              ,value  :: generic_id
    type(FT_Epetra_RowMatrix_ID_t),pointer:: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraRowMatrix = local_ptr
  end function
 
  subroutine RowMatrix_assign_ID(lhs,rhs)
    use ForTrilinos_enums
    type(FT_Epetra_RowMatrix_ID_t) ,intent(in)   :: rhs
    class(Epetra_RowMatrix)        ,intent(inout):: lhs
    lhs%RowMatrix_id=rhs
  end subroutine
 
  subroutine invalidate_EpetraRowMatrix_ID(this)
    class(Epetra_RowMatrix) ,intent(inout) :: this
    this%RowMatrix_id%table = FT_Invalid_ID
    this%RowMatrix_id%index = FT_Invalid_Index 
    this%RowMatrix_id%is_const = FT_FALSE
  end subroutine
 
  subroutine ctrilinos_delete_EpetraRowMatrix(this)
    class(Epetra_RowMatrix) ,intent(inout) :: this
    call Epetra_RowMatrix_Destroy( this%RowMatrix_id )
  end subroutine
end module 
