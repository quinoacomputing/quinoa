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


module FAmesos
  use ForTrilinos_enums ,only: FT_Amesos_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use iso_c_binding     ,only: c_int,c_double,c_char
  use foramesos
  implicit none
  private          ! Hide everything by default
  public :: Amesos ! Expose type/constructors/methods

  type ,extends(universal) :: Amesos 
    private
    type(FT_Amesos_ID_t) :: Amesos_id 
  contains
     ! Constructors
     procedure ,private :: from_struct_
     procedure ,private :: create_
     generic :: Amesos_ => from_struct_,create_
     ! Developers only
     procedure         :: invalidate_id => invalidate_Amesos_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_Amesos
     procedure         :: get_Amesos_ID 
     procedure ,nopass :: alias_Amesos_ID
     procedure         :: generalize 
  end type

   interface Amesos ! constructors
     module procedure from_struct,create
   end interface

contains

  subroutine from_struct_(this,id)
    class(Amesos) ,intent(out) :: this
    type(FT_Amesos_ID_t) ,intent(in) :: id
    this%Amesos_id = id
    call this%register_self
  end subroutine
  
  function from_struct(id) result(new_Amesos)
    type(Amesos) :: new_Amesos
    type(FT_Amesos_ID_t) ,intent(in) :: id
    new_Amesos%Amesos_id = id
    call new_Amesos%register_self
  end function
  
  subroutine create_(this)
    class(Amesos) ,intent(out) :: this
    call this%Amesos_(Amesos_Create())
  end subroutine

  function create() result(new_Amesos)
    type(Amesos) :: new_Amesos
    call new_Amesos%Amesos_()
  end function

  type(FT_Amesos_ID_t) function get_Amesos_ID(this)
    class(Amesos) ,intent(in) :: this 
    get_Amesos_ID=this%Amesos_id
  end function
  
  type(FT_Amesos_ID_t) function alias_Amesos_ID(generic_id)
    use ForTrilinos_table_man,only: CT_Alias
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t,FT_Amesos_ID
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Amesos_ID),stat=status)
    ierr=error(status,'FAmesos:alias_Amesos_ID')
    call ierr%check_success()
    alias_Amesos_ID=degeneralize_Amesos(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only : c_loc
   class(Amesos) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%Amesos_id))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Amesos) ,intent(in) ,target :: this
   ! generalize = Amesos_Generalize ( this%Amesos_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Amesos_ID_t) function degeneralize_Amesos(generic_id) bind(C)
  ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Amesos_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                      ,value   :: generic_id
    type(FT_Amesos_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_Amesos = local_ptr
  ! ____ Use for ForTrilinos function implementation ______

  ! ____ Use for CTrilinos function implementation ______
  !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
  !degeneralize_Amesos = AztecOO_Degeneralize(generic_id)
  ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine invalidate_Amesos_ID(this)
    class(Amesos),intent(inout) :: this
    this%Amesos_id%table = FT_Invalid_ID
    this%Amesos_id%index = FT_Invalid_Index 
    this%Amesos_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_Amesos(this)
    class(Amesos),intent(inout) :: this
    call Amesos_Destroy( this%Amesos_id ) 
  end subroutine

end module 

