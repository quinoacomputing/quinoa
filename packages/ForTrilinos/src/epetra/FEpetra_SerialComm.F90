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

#include "ForTrilinos_config.h"
module FEpetra_SerialComm
  use ForTrilinos_enums ,only : FT_Epetra_Comm_ID,FT_Epetra_SerialComm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_error ,only : error
  use ForTrilinos_assertion_utility
  use FEpetra_Comm      ,only : Epetra_Comm
  use iso_c_binding     ,only : c_int,c_long,c_double,c_char
  use forepetra
  implicit none
  private                     ! Hide everything by default
  public :: Epetra_SerialComm ! Expose type/constructors/methods

  type ,extends(Epetra_Comm)        :: Epetra_SerialComm 
    private
    type(FT_Epetra_SerialComm_ID_t) :: SerialComm_id 
  contains
     !Constructors
     procedure ,private :: duplicate_
     procedure ,private :: create_
     procedure ,private :: from_comm_ 
     procedure ,private :: from_struct_
     generic :: Epetra_SerialComm_ => create_,from_comm_,from_struct_,duplicate_
     !Barrier Methods
     procedure         :: barrier
     !Broadcast Methods
     procedure :: broadcast_double
     procedure :: broadcast_int
     procedure :: broadcast_long
     procedure :: broadcast_char
     !Gather Methods
     procedure :: gather_double
     procedure :: gather_int
     procedure :: gather_long
     !Sum Methods
     procedure :: sum_double
     procedure :: sum_int
     procedure :: sum_long
     !Max/Min Methods
     procedure :: max_double
     procedure :: max_int
     procedure :: max_long
     procedure :: min_double
     procedure :: min_int
     procedure :: min_long
     !Parallel Prefix Methods
     procedure :: ScanSum_double
     procedure :: ScanSum_int
     procedure :: ScanSum_long
     !Attribute Accessor Methods
     procedure :: MyPID
     procedure :: NumProc
     !Developers only  -- to be called by developers from other ForTrilinos modules, not by end applications:
     procedure          :: invalidate_id => invalidate_EpetraSerialComm_ID 
     procedure          :: ctrilinos_delete => ctrilinos_delete_EpetraSerialComm
     procedure          :: get_EpetraSerialComm_ID 
     procedure ,nopass  :: alias_EpetraSerialComm_ID
     procedure          :: generalize 
  end type

   interface Epetra_SerialComm ! constructors
    !User interface -- constructors for use by end applications:
    module procedure create,duplicate
    !Developers only -- to be called by developers from other ForTrilinos modules, not by end applications:
    module procedure from_struct, from_comm
   end interface

contains
  ! Common type-bound procedures begin here: all ForTrilinos child classes of the API have procedures analogous to these.
  ! The function from_struct must construct a struct id for the extended parent class
  ! The function from_struct should be called by developers only from other ForTrilinos modules, never by end applications:

  subroutine from_struct_(this,id)
    class(Epetra_SerialComm) ,intent(out) :: this
    type(FT_Epetra_SerialComm_ID_t) ,intent(in) :: id
    this%SerialComm_id = id
    call this%set_EpetraComm_ID(this%alias_EpetraComm_ID(this%generalize()))
    call this%register_self
  end subroutine

  function from_struct(id) result(new_Epetra_SerialComm)
    type(FT_Epetra_SerialComm_ID_t) ,intent(in) :: id
    type(Epetra_SerialComm) :: new_Epetra_SerialComm
    call new_Epetra_SerialComm%Epetra_SerialComm_(id)
  end function

  subroutine from_comm_(this,id)
    class(Epetra_SerialComm) ,intent(out) :: this 
    type(FT_Epetra_Comm_ID_t) ,intent(in) :: id
    call this%set_EpetraComm_ID(id)
    this%SerialComm_id = this%alias_EpetraSerialComm_ID(this%generalize_EpetraComm())
    call this%register_self
  end subroutine

  function from_comm(id) result(new_Epetra_SerialComm)
    type(Epetra_SerialComm) :: new_Epetra_SerialComm
    type(FT_Epetra_Comm_ID_t) ,intent(in) :: id
    call new_Epetra_SerialComm%Epetra_SerialComm_(id)
  end function

  subroutine create_(this)
    class(Epetra_SerialComm) ,intent(out) :: this 
    call this%Epetra_SerialComm_(Epetra_SerialComm_Create())
  end subroutine

  function create() result(new_Epetra_SerialComm)
    type(Epetra_SerialComm) :: new_Epetra_SerialComm
    call new_Epetra_SerialComm%Epetra_SerialComm_()
  end function

  subroutine duplicate_(this,copy)
    class(Epetra_SerialComm) ,intent(in) :: this
    type(Epetra_SerialComm) ,intent(out) :: copy
    call copy%Epetra_SerialComm_(Epetra_SerialComm_Duplicate(this%SerialComm_id))
  end subroutine

  function duplicate(original) result (new_Epetra_SerialComm)
    type(Epetra_SerialComm) ,intent(in) :: original
    type(Epetra_SerialComm) :: new_Epetra_SerialComm
    call original%Epetra_SerialComm_(new_Epetra_SerialComm)
  end function

  !----------------- Struct access ---------------------------------------------

  type(FT_Epetra_SerialComm_ID_t) function get_EpetraSerialComm_ID(this)
    class(Epetra_SerialComm) ,intent(in) :: this 
    get_EpetraSerialComm_ID=this%SerialComm_id
  end function
  
 !----------------- Type casting ---------------------------------------------

  type(FT_Epetra_SerialComm_ID_t) function alias_EpetraSerialComm_ID(generic_id)
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: FT_Epetra_SerialComm_ID,ForTrilinos_Universal_ID_t
    use ForTrilinos_table_man,only: CT_Alias
    type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(Fortrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_SerialComm_ID),stat=status)
    ierr=error(status,'FEpetra_SerialComm:alias_EpetraSerialComm_ID')
    call ierr%check_success()
    alias_EpetraSerialComm_ID=degeneralize_EpetraSerialComm(c_loc(alias_id))
  end function


  type(ForTrilinos_Universal_ID_t) function generalize(this)
    use ForTrilinos_utils ,only: generalize_all
    use iso_c_binding ,only : c_loc
    class(Epetra_SerialComm) ,intent(in) ,target :: this
    generalize = generalize_all( c_loc(this%SerialComm_id) )
  end function
 
 type(FT_Epetra_SerialComm_ID_t) function degeneralize_EpetraSerialComm(generic_id) bind(C)
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_SerialComm_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                     ,value   :: generic_id
    type(FT_Epetra_SerialComm_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraSerialComm = local_ptr
  end function
 
  subroutine barrier(this)
    class(Epetra_SerialComm) ,intent(in) :: this
    call Epetra_SerialComm_Barrier(this%SerialComm_id)
  end subroutine
 
  subroutine broadcast_double(this,MyVals,root,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(inout) :: MyVals
    integer(c_int)               ,intent(in)    :: root
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_Broadcast_Double(this%SerialComm_id,MyVals,size(MyVals),root)
    if (present(err)) err=error(error_out,'Epetra_SerialComm%broadcast_double: failed.')
  end subroutine

  subroutine broadcast_int(this,MyVals,root,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(inout) :: MyVals
    integer(c_int)               ,intent(in)    :: root
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_Broadcast_Int(this%SerialComm_id,MyVals,size(MyVals),root)
    if (present(err)) err=error(error_out,'Epetra_SerialComm%broadcast_int: failed.')
  end subroutine

  subroutine broadcast_long(this,MyVals,root,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long),dimension(:) ,intent(inout) :: MyVals
    integer(c_int)               ,intent(in)    :: root
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_Broadcast_Long(this%SerialComm_id,MyVals,size(MyVals),root)
    if (present(err)) err=error(error_out,'Epetra_SerialComm%broadcast_long: failed.')
  end subroutine
 
  subroutine broadcast_char(this,MyVals,root,err)
    class(Epetra_SerialComm)           ,intent(in)    :: this
    character(kind=c_char),dimension(:),intent(inout) :: MyVals
    integer(c_int)                     ,intent(in)    :: root
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_Broadcast_Char(this%SerialComm_id,MyVals,size(MyVals),root)
    if (present(err)) err=error(error_out,'Epetra_SerialComm%broadcast_char: failed.')
  end subroutine
  
 subroutine gather_double(this,MyVals,AllVals,err)
   class(Epetra_SerialComm)     ,intent(in)    :: this
   real(c_double), dimension(:) ,intent(in)    :: MyVals
   real(c_double), dimension(:) ,intent(inout) :: AllVals
   type(error) ,optional, intent(inout) :: err
   integer(c_int)     :: error_out
   call assert(size(AllVals)==this%NumProc()*size(MyVals)&
     ,error_message('GatherAll: AllVals must be of size NumProc*size(MyVals)' ))
   error_out = Epetra_SerialComm_GatherAll_Double(this%SerialComm_id,MyVals,AllVals,size(MyVals))
   if (present(err)) err=error(error_out,'Epetra_SerialComm%gather_double: failed.')
  end subroutine

  subroutine gather_int(this,MyVals,AllVals,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: MyVals
    integer(c_int), dimension(:) ,intent(inout) :: AllVals
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    call assert(size(AllVals)==this%NumProc()*size(MyVals)&
      ,error_message('GatherAll: AllVals must be of size NumProc*size(MyVals)'))
    error_out = Epetra_SerialComm_GatherAll_Int(this%SerialComm_id,MyVals,AllVals,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%gather_int: failed.')
  end subroutine

  subroutine gather_long(this,MyVals,AllVals,err)
    class(Epetra_SerialComm)      ,intent(in)    :: this
    integer(c_long), dimension(:) ,intent(in)    :: MyVals
    integer(c_long), dimension(:) ,intent(inout) :: AllVals
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    call assert(size(AllVals)==this%NumProc()*size(MyVals)&
      ,error_message('gather_long: AllVals must be of size NumProc*size(MyVals)' ))
    error_out = Epetra_SerialComm_GatherAll_Long(this%SerialComm_id,MyVals,AllVals,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%gather_long: failed.')
  end subroutine

  subroutine sum_double(this,PartialSums,GlobalSums,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialSums
    real(c_double), dimension(:) ,intent(inout) :: GlobalSums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_SumAll_Double(this%SerialComm_id,PartialSums,GlobalSums,size(PartialSums))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%sum_double: failed.')
  end subroutine

  subroutine sum_int(this,PartialSums,GlobalSums,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialSums
    integer(c_int), dimension(:) ,intent(inout) :: GlobalSums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_SumAll_Int(this%SerialComm_id,PartialSums,GlobalSums,size(PartialSums))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%sum_int: failed.')
  end subroutine

  subroutine sum_long(this,PartialSums,GlobalSums,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialSums
     integer(c_long), dimension(:),intent(inout) :: GlobalSums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_SumAll_Long(this%SerialComm_id,PartialSums,GlobalSums,size(partialSums))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%sum_long: failed.')
  end subroutine
  
  subroutine max_double(this,PartialMaxs,GlobalMaxs,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialMaxs
    real(c_double), dimension(:) ,intent(inout) :: GlobalMaxs
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_MaxAll_Double(this%SerialComm_id,PartialMaxs,GlobalMaxs,size(PartialMaxs))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%max_double: failed.')
  end subroutine

  subroutine max_int(this,PartialMaxs,GlobalMaxs,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialMaxs
    integer(c_int), dimension(:) ,intent(inout) :: GlobalMaxs
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_MaxAll_Int(this%SerialComm_id,PartialMaxs,GlobalMaxs,size(PartialMaxs))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%max_int: failed.')
  end subroutine

  subroutine max_long(this,PartialMaxs,GlobalMaxs,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialMaxs
    integer(c_long), dimension(:),intent(inout) :: GlobalMaxs
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_MaxAll_Long(this%SerialComm_id,PartialMaxs,GlobalMaxs,size(PartialMaxs))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%max_long: failed.')
  end subroutine
  
  subroutine min_double(this,PartialMins,GlobalMins,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialMins
    real(c_double), dimension(:) ,intent(inout) :: GlobalMins
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_MinAll_Double(this%SerialComm_id,PartialMins,GlobalMins,size(PartialMins))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%min_double: failed.')
  end subroutine

  subroutine min_int(this,PartialMins,GlobalMins,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialMins
    integer(c_int), dimension(:) ,intent(inout) :: GlobalMins
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_MinAll_Int(this%SerialComm_id,PartialMins,GlobalMins,size(PartialMins))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%min_int: failed.')
  end subroutine

  subroutine min_long(this,PartialMins,GlobalMins,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialMins
    integer(c_long), dimension(:),intent(inout) :: GlobalMins
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_MinAll_Long(this%SerialComm_id,PartialMins,GlobalMins,size(PartialMins))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%min_long: failed.')
  end subroutine

  subroutine ScanSum_double(this,MyVals,scan_sums,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: MyVals 
    real(c_double), dimension(:) ,intent(inout) :: scan_sums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_ScanSum_Double(this%SerialComm_id,MyVals,scan_sums,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%ScanSum_double: failed.')
  end subroutine

  subroutine ScanSum_int(this,MyVals,scan_sums,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: MyVals 
    integer(c_int), dimension(:) ,intent(inout) :: scan_sums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_ScanSum_Int(this%SerialComm_id,MyVals,scan_sums,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%ScanSum_int: failed.')
  end subroutine

  subroutine ScanSum_long(this,MyVals,scan_sums,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: MyVals 
     integer(c_long), dimension(:),intent(inout) :: scan_sums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_SerialComm_ScanSum_Long(this%SerialComm_id,MyVals,scan_sums,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_SerialComm%ScanSum_long: failed.')
  end subroutine

  integer(c_int) function MyPID(this)
    class(Epetra_SerialComm)     , intent(in) :: this
    MyPID=Epetra_SerialComm_MyPID(this%SerialComm_id)
  end function

  integer(c_int) function NumProc(this)
    class(Epetra_SerialComm)     , intent(in) :: this
    NumProc=Epetra_SerialComm_NumProc(this%SerialComm_id)
  end function

  !__________ Garbage collection __________________________________________________ 

  subroutine invalidate_EpetraSerialComm_ID(this)
    class(Epetra_SerialComm),intent(inout) :: this
    call this%invalidate_EpetraComm_ID
    this%SerialComm_id%table=FT_Invalid_ID
    this%SerialComm_id%index=FT_Invalid_Index
    this%SerialComm_id%is_const=FT_FALSE
  end subroutine
  
  subroutine ctrilinos_delete_EpetraSerialComm(this)
    class(Epetra_SerialComm) ,intent(inout) :: this
    call this%ctrilinos_delete_EpetraComm()
    call Epetra_SerialComm_Destroy(this%SerialComm_id)
  end subroutine

end module 
