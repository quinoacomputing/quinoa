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


module FEpetra_Map
  use ForTrilinos_enums ,only: FT_Epetra_BlockMap_ID,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_error
  use FEpetra_Comm       ,only: Epetra_Comm
  use FEpetra_BlockMap   ,only: Epetra_BlockMap
  use iso_c_binding      ,only: c_int
  use forepetra
  implicit none
  private                     ! Hide everything by default
  public :: Epetra_Map        ! Expose type/constructors/methods

  type, extends(Epetra_BlockMap) :: Epetra_Map 
    private
    type(FT_Epetra_Map_ID_t) :: map_id
  contains
     !Constructors
     procedure ,private :: create__
     procedure ,private :: create_linear__
     procedure ,private :: duplicate__
     procedure ,private :: create_arbitrary__
     !Developers only  -- to be called by developers from other ForTrilinos modules, not by end applications:
     procedure ,private :: from_struct__
     generic :: Epetra_Map_ => create__,from_struct__,create_linear__,duplicate__,create_arbitrary__
     procedure         :: invalidate_id => invalidate_EpetraMap_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraMap
     procedure         :: get_EpetraMap_ID 
     procedure ,nopass :: alias_EpetraMap_ID
     procedure         :: generalize 
  end type

   interface Epetra_Map ! constructors
     !User interface -- constructors for use by end applications:
     module procedure create,duplicate,create_linear,create_arbitrary
     !Developers only -- to be called by developers from other ForTrilinos modules, not by end applications:
     module procedure from_struct
   end interface
 
contains
  ! Common type-bound procedures begin here: all ForTrilinos child classes of the API have procedures analogous to these.
  ! The subroutine from_struct_ must construct a struct id for the extended parent class.  The procedures 
  ! from_struct_ and from_struct should be called by developers only from other ForTrilinos modules, never by end applications:

  subroutine from_struct__(this,id)
    class(Epetra_Map) ,intent(out) :: this
    type(FT_Epetra_Map_ID_t),intent(in) :: id
    this%map_id = id
    this%Epetra_BlockMap=Epetra_BlockMap(this%alias_EpetraBlockMap_ID(this%generalize()))
    call this%register_self()
  end subroutine

  function from_struct(id) result(new_Epetra_Map)
    type(Epetra_Map) :: new_Epetra_Map
    type(FT_Epetra_Map_ID_t),intent(in) :: id
    call new_Epetra_Map%Epetra_Map_(id)
  end function

  ! All additional constructors should take two steps: (1) obtain a struct ID by invoking a procedural binding and then (2) pass
  ! this ID to from_struct to initialize the constructed object's ID component and register the object for reference counting.

  subroutine create__(this,Num_GlobalElements,IndexBase,comm)
    class(Epetra_Map) ,intent(out) :: this
    integer(c_int) ,intent(in) :: Num_GlobalElements,IndexBase
    class(Epetra_Comm) ,intent(in) :: comm
    call this%Epetra_Map_(Epetra_Map_Create(Num_GlobalElements,IndexBase,comm%get_EpetraComm_ID()))
  end subroutine

  function create(Num_GlobalElements,IndexBase,comm) result(new_Epetra_Map)
    type(Epetra_Map) :: new_Epetra_Map
    integer(c_int) ,intent(in) :: Num_GlobalElements,IndexBase
    class(Epetra_Comm) ,intent(in) :: comm
    call new_Epetra_Map%Epetra_Map_(Num_GlobalElements,IndexBase,comm)
  end function

  subroutine create_linear__(this,Num_GlobalElements,Num_MyElements,IndexBase,comm)
    class(Epetra_Map) ,intent(out) :: this
    integer(c_int) ,intent(in) :: Num_GlobalElements,Num_MyElements,IndexBase
    class(Epetra_Comm) ,intent(in) :: comm
    call this%Epetra_Map_(Epetra_Map_Create_Linear(Num_GlobalElements,Num_MyElements,IndexBase,comm%get_EpetraComm_ID()))
  end subroutine
  
  function create_linear(Num_GlobalElements,Num_MyElements,IndexBase,comm) result(new_Epetra_Map)
    type(Epetra_Map) :: new_Epetra_Map
    integer(c_int) ,intent(in) :: Num_GlobalElements,Num_MyElements,IndexBase
    class(Epetra_Comm) ,intent(in) :: comm
    call new_Epetra_Map%Epetra_Map_(Num_GlobalElements,Num_MyElements,IndexBase,comm)
  end function
  
  subroutine create_arbitrary__(this,Num_GlobalElements,My_GlobalElements,IndexBase,comm)
    class(Epetra_Map) ,intent(out) :: this
    integer(c_int) ,intent(in)              :: Num_GlobalElements,IndexBase
    integer(c_int) ,intent(in) ,dimension(:),allocatable:: My_GlobalElements
    class(Epetra_Comm) ,intent(in) :: comm
    call this% & 
    Epetra_Map_(Epetra_Map_Create_Arbitrary(Num_GlobalElements,size(My_GlobalElements),&
               My_GlobalElements,IndexBase,comm%get_EpetraComm_ID()))
  end subroutine
 
  function create_arbitrary(Num_GlobalElements,My_GlobalElements,IndexBase,comm) result(new_Epetra_Map)
    type(Epetra_Map) :: new_Epetra_Map
    integer(c_int) ,intent(in)              :: Num_GlobalElements,IndexBase
    integer(c_int) ,intent(in) ,dimension(:),allocatable:: My_GlobalElements
    class(Epetra_Comm) ,intent(in) :: comm
    call new_Epetra_Map%Epetra_Map_(Num_GlobalElements,My_GlobalElements,IndexBase,comm)
  end function
 
  subroutine duplicate__(this,copy)
    class(Epetra_Map) ,intent(in) :: this
    type(Epetra_Map) ,intent(out) :: copy
    call copy%Epetra_Map_(Epetra_Map_Duplicate(this%map_id))
  end subroutine

  type(Epetra_Map) function duplicate(original)
    type(Epetra_Map) ,intent(in) :: original
    call original%Epetra_Map_(duplicate)
  end function

  !----------------- Struct access ---------------------------------------------

  type(FT_Epetra_Map_ID_t) function get_EpetraMap_ID(this)
    class(Epetra_Map) ,intent(in) :: this 
    get_EpetraMap_ID=this%map_id
  end function

  !----------------- Type casting ---------------------------------------------
 
  type(FT_Epetra_Map_ID_t) function alias_EpetraMap_ID(generic_id)
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t, FT_Epetra_Map_ID
    use ForTrilinos_table_man,only: CT_Alias
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Map_ID),stat=status)
    ierr=error(status,'FEpetra_Map:alias_EpetraMap_ID')
    call ierr%check_success()
    alias_EpetraMap_ID=degeneralize_EpetraMap(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   use ForTrilinos_utils ,only: generalize_all 
   use iso_c_binding     ,only: c_loc
   class(Epetra_Map) ,intent(in) ,target :: this
   generalize =generalize_all(c_loc(this%map_id))
  end function
  
  type(FT_Epetra_Map_ID_t) function degeneralize_EpetraMap(generic_id) bind(C)
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_Map_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)              ,value   :: generic_id
    type(FT_Epetra_Map_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraMap = local_ptr
    local_ptr => null()
  end function

 !__________ Garbage collection __________________________________________________

  subroutine invalidate_EpetraMap_ID(this)
    class(Epetra_Map),intent(inout) :: this
    call this%Epetra_BlockMap%invalidate_id
    this%Map_id%table = FT_Invalid_ID
    this%Map_id%index = FT_Invalid_Index 
    this%Map_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraMap(this)
    class(Epetra_Map),intent(inout) :: this
    call Epetra_Map_Destroy( this%map_id ) 
  end subroutine

end module 

