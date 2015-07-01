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


module FEpetra_SrcDistObject
  use ForTrilinos_enums !,only: FT_Epetra_SrcDistObject_ID_t,ForTrilinos_Universal_ID_t,FT_Epetra_BlockMap_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only : universal
  use ForTrilinos_error
  use FEpetra_BlockMap, only: Epetra_BlockMap
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: Epetra_SrcDistObject ! Expose type/methods

  type ,abstract ,extends(universal) :: Epetra_SrcDistObject
   private
   type(FT_Epetra_SrcDistObject_ID_t)         :: SrcDistObject_id 
  contains
    ! Developers only
    procedure                          :: invalidate_EpetraSrcDistObject_ID
    procedure                          :: ctrilinos_delete_EpetraSrcDistObject
    procedure                          :: get_EpetraSrcDistObject_ID
    procedure                          :: set_EpetraSrcDistObject_ID
    procedure                 ,nopass  :: alias_EpetraSrcDistObject_ID
    procedure ,non_overridable         :: generalize_EpetraSrcDistObject
    !procedure(map_interface),deferred  :: map
  end type
  
  
  !abstract interface
  ! We need a CTrilinos function that tells us the type of the underlying object
  !function map_interface(this) result(map)  
  !  import :: Epetra_BlockMap, Epetra_SrcDistObject
  !  class(Epetra_BlockMap), allocatable :: map  
  !  class(Epetra_SrcDistObject), intent(in) :: this
  !end function
  !end interface
 
  contains
  type(FT_Epetra_SrcDistObject_ID_t) function get_EpetraSrcDistObject_ID(this)
    class(Epetra_SrcDistObject) ,intent(in) :: this
    get_EpetraSrcDistObject_ID = this%SrcDistObject_id
  end function
  
  subroutine set_EpetraSrcDistObject_ID(this,id)
    class(Epetra_SrcDistObject)        ,intent(inout) :: this
    type(FT_Epetra_SrcDistObject_ID_t) ,intent(in)    :: id 
    this%SrcDistObject_id=id
  end subroutine 
  
  type(FT_Epetra_SrcDistObject_ID_t) function alias_EpetraSrcDistObject_ID(generic_id)
    use iso_c_binding, only : c_loc,c_int
    use ForTrilinos_table_man
    use ForTrilinos_enums
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_SrcDistObject_ID),stat=status)
    ierr=error(status,'FEpetra_SrcDistObject:alias_EpetraSrcDistObject_ID')
    call ierr%check_success()
    alias_EpetraSrcDistObject_ID=degeneralize_EpetraSrcDistObject(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize_EpetraSrcDistObject(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_SrcDistObject) ,intent(in) ,target :: this
   generalize_EpetraSrcDistObject = generalize_all( c_loc(this%SrcDistObject_id) )
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_SrcDistObject) ,intent(in) ,target :: this
   ! generalize_EpetraSrcDistObject = Epetra_RowMatrix_Generalize ( this%SrcDistObject_id )
   ! ____ Use for CTrilinos function implementation ______
  end function
  
  type(FT_Epetra_SrcDistObject_ID_t) function degeneralize_EpetraSrcDistObject(generic_id) bind(C)
    !use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_SrcDistObject_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)              ,value  :: generic_id
    type(FT_Epetra_SrcDistObject_ID_t),pointer:: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraSrcDistObject = local_ptr
  end function
 
  subroutine SrcDistObject_assign_ID(lhs,rhs)
    use ForTrilinos_enums
    type(FT_Epetra_SrcDistObject_ID_t) ,intent(in)   :: rhs
    class(Epetra_SrcDistObject)        ,intent(inout):: lhs
    lhs%SrcDistObject_id=rhs
  end subroutine
 
  subroutine invalidate_EpetraSrcDistObject_ID(this)
    class(Epetra_SrcDistObject) ,intent(inout) :: this
    this%SrcDistObject_id%table = FT_Invalid_ID
    this%SrcDistObject_id%index = FT_Invalid_Index 
    this%SrcDistObject_id%is_const = FT_FALSE
  end subroutine
 
  subroutine ctrilinos_delete_EpetraSrcDistObject(this)
    class(Epetra_SrcDistObject) ,intent(inout) :: this
    call Epetra_SrcDistObject_Destroy( this%SrcDistObject_id )
  end subroutine
end module 
