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
module FEpetra_MpiComm
#ifdef HAVE_MPI
  use ForTrilinos_enums ,only: FT_Epetra_MpiComm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_error ,only :error
  use ForTrilinos_external_utils
  use ForTrilinos_assertion_utility
  use FEpetra_Comm      ,only: Epetra_Comm
  use iso_c_binding     ,only: c_int,c_double,c_long,c_char
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: Epetra_MpiComm ! Expose type/methods

  type ,extends(Epetra_Comm) :: Epetra_MpiComm
    private
    type(FT_Epetra_MpiComm_ID_t) :: MpiComm_id  
  contains
   ! Constructors
    procedure ,private :: from_scratch_
    procedure ,private :: duplicate_
    procedure ,private :: from_struct_
    generic :: Epetra_MpiComm_ => from_scratch_,duplicate_,from_struct_
    !Barrier Method
    procedure         :: barrier
    !Broadcast Method
    procedure         :: broadcast_double
    procedure         :: broadcast_int
    procedure         :: broadcast_long
    procedure         :: broadcast_char
    !Gather Methods
    procedure         :: gather_double
    procedure         :: gather_int
    procedure         :: gather_long
    !Sum Methods
    procedure         :: sum_double
    procedure         :: sum_int
    procedure         :: sum_long
    !Max/Min Methods
    procedure         :: max_double
    procedure         :: max_int
    procedure         :: max_long
    procedure         :: min_double
    procedure         :: min_int
    procedure         :: min_long
    !Parallel Prefix Methods
    procedure         :: ScanSum_double
    procedure         :: ScanSum_int
    procedure         :: ScanSum_long
    !Attribute Accessor Methods
    procedure         :: MyPID
    procedure         :: NumProc
    !The following procedures are for ForTrilinos developers only (not for users):
    procedure         :: invalidate_id => invalidate_EpetraMpiComm_ID
    procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraMpiComm
    procedure         :: get_EpetraMpiComm_ID
    procedure ,nopass :: alias_EpetraMpiComm_ID
    procedure         :: generalize
  end type
  
  interface Epetra_MpiComm ! constructor functions
    module procedure from_scratch,duplicate
    ! The following procedures are for ForTrilinos developers only (not for users):
    module procedure from_struct,from_comm
  end interface

contains

  ! Every constructor comes in two flavors: a function and a subroutine.  Each function is a simple pass-through wrapper
  ! for the corresponding subroutine.  This structure serves a useful purpose for unit testing: invoking the function tests
  ! both flavors.  Although some situations necessitate the function (e.g., passing a function invocation as an argument to
  ! a procedure), the subroutine is preferred because it incurs considerably less overhead related to the ForTrilinos 
  ! reference-counting architecture. 

  ! Every constructor must register the constructed object for reference-counting by invoking register_self on the object  
  ! or by passing the object to a procedure that invokes register_self on the object.

  subroutine from_struct_(this,id)
    class(Epetra_MpiComm) ,intent(out) :: this
    type(FT_Epetra_MpiComm_ID_t) ,intent(in) :: id
    this%MpiComm_id = id
    call this%set_EpetraComm_ID(this%alias_EpetraComm_ID(this%generalize()))
    call this%register_self
  end subroutine

  function from_struct(id) result(new_Epetra_MpiComm)
    type(Epetra_MpiComm) :: new_Epetra_MpiComm 
    type(FT_Epetra_MpiComm_ID_t) ,intent(in) :: id
    call new_Epetra_MpiComm%Epetra_MpiComm_(id)
  end function
 
 type(Epetra_MpiComm) function from_comm(id)
   type(FT_Epetra_Comm_ID_t) ,intent(in) :: id
   call from_comm%set_EpetraComm_ID(id)
   from_comm%MpiComm_id = from_comm%alias_EpetraMpiComm_ID(from_comm%generalize_EpetraComm())
   call from_comm%register_self
  end function

  ! Original C++ prototype:
  ! Epetra_MpiComm(MPI_Comm comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm );

  subroutine from_scratch_(this,comm)
   class(Epetra_MpiComm) ,intent(out) :: this
   integer(c_int) ,intent(in) :: comm
   call this%Epetra_MpiComm_(Epetra_MpiComm_Fortran_Create(comm))
  end subroutine

  function from_scratch(comm) result(new_Epetra_MpiComm)
    type(Epetra_MpiComm) :: new_Epetra_MpiComm
    integer(c_int) ,intent(in) :: comm
    call new_Epetra_MpiComm%Epetra_MpiComm_(comm) 
  end function

  ! Original C++ prototype:
  ! Epetra_MpiComm(const Epetra_MpiComm & Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( CT_Epetra_MpiComm_ID_t CommID );

  subroutine duplicate_(this,copy)
    class(Epetra_MpiComm) ,intent(in) :: this
    type(Epetra_MpiComm) ,intent(out) :: copy
    call copy%Epetra_MpiComm_(Epetra_MpiComm_Duplicate(this%MpiComm_id))
  end subroutine

  type(Epetra_MpiComm) function duplicate(original)
    type(Epetra_MpiComm) ,intent(in) :: original
    call original%Epetra_MpiComm_(duplicate) 
  end function

  type(FT_Epetra_MpiComm_ID_t) function get_EpetraMpiComm_ID(this)
    class(Epetra_MpiComm) ,intent(in) :: this
    get_EpetraMpiComm_ID=this%MpiComm_id
  end function
 
  type(FT_Epetra_MpiComm_ID_t) function alias_EpetraMpiComm_ID(generic_id)
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: FT_Epetra_MpiComm_ID,ForTrilinos_Universal_ID_t
    use ForTrilinos_table_man,only: CT_Alias
    type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(Fortrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_MpiComm_ID),stat=status)
    ierr=error(status,'FEpetra_MpiComm:alias_EpetraMpiComm_ID')
    call ierr%check_success()
    alias_EpetraMpiComm_ID=degeneralize_EpetraMpiComm(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_MpiComm) ,intent(in) ,target :: this
   generalize = generalize_all( c_loc(this%MpiComm_id) )
   ! ____ Use for ForTrilinos function implementation ______
  
   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_MpiComm) ,intent(in) ,target :: this
   ! generalize = Epetra_MpiComm_Generalize ( this%MpiComm_id )
   ! ____ Use for CTrilinos function implementation ______
  end function

  type(FT_Epetra_MpiComm_ID_t) function degeneralize_EpetraMpiComm(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_MpiComm_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                     ,value   :: generic_id
    type(FT_Epetra_MpiComm_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraMpiComm = local_ptr
   ! ____ Use for ForTrilinos function implementation ______
  
   ! ____ Use for CTrilinos function implementation ______
   ! type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
   ! degeneralize_EpetraMpiComm = Epetra_MpiComm_Degeneralize( generic_id )
   ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine barrier(this)
   class(Epetra_MpiComm) ,intent(in) :: this
   call Epetra_MpiComm_Barrier(this%MpiComm_id)
  end subroutine

  subroutine broadcast_double(this,MyVals,root,err)
    class(Epetra_MpiComm)       ,intent(in)    :: this
    real(c_double) ,dimension(:),intent(inout) :: MyVals
    integer(c_int)              ,intent(in)    :: root
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out=Epetra_MpiComm_Broadcast_Double(this%MpiComm_id,MyVals,size(MyVals),root)
    if (present(err)) err=error(error_out,'Epetra_MpiComm%broadcast_double: failed.')
  end subroutine

  subroutine broadcast_int(this,MyVals,root,err)
    class(Epetra_MpiComm)       ,intent(in)    :: this
    integer(c_int) ,dimension(:),intent(inout) :: MyVals
    integer(c_int)              ,intent(in)    :: root
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out=Epetra_MpiComm_Broadcast_Int(this%MpiComm_id,MyVals,size(MyVals),root)
    if (present(err)) err=error(error_out,'Epetra_MpiComm%broadcast_int: failed.')
  end subroutine

  subroutine broadcast_long(this,MyVals,root,err)
    class(Epetra_MpiComm)       ,intent(in)    :: this
    integer(c_long),dimension(:),intent(inout) :: MyVals
    integer(c_int)              ,intent(in)    :: root
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out=Epetra_MpiComm_Broadcast_Long(this%MpiComm_id,MyVals,size(MyVals),root)
    if (present(err)) err=error(error_out,'Epetra_MpiComm%broadcast_long: failed.')
  end subroutine

  subroutine broadcast_char(this,MyVals,root,err)
    class(Epetra_MpiComm)              ,intent(in)    :: this
    character(kind=c_char),dimension(:),intent(inout) :: MyVals
    integer(c_int)                     ,intent(in)    :: root
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out=Epetra_MpiComm_Broadcast_Char(this%MpiComm_id,MyVals,size(MyVals),root)
    if (present(err)) err=error(error_out,'Epetra_MpiComm%broadcast_char: failed.')
  end subroutine

  subroutine gather_double(this,MyVals,AllVals,err)
    class(Epetra_MpiComm)      ,intent(in)    :: this
    real(c_double),dimension(:),intent(in)    :: MyVals
    real(c_double),dimension(:),intent(inout) :: AllVals
    type(error) ,optional      ,intent(inout) :: err
    integer(c_int)     :: error_out
    call assert(size(AllVals)==this%NumProc()*size(MyVals)&
        ,error_message('GatherAll: AllVals must be of size NumProc*size(MyVals)'))
    error_out = Epetra_MpiComm_GatherAll_Double(this%MpiComm_id,MyVals,AllVals,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%gather_double: failed.')
  end subroutine

  subroutine gather_int(this,MyVals,AllVals,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: MyVals
    integer(c_int), dimension(:) ,intent(inout) :: AllVals
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    call assert(size(AllVals)==this%NumProc()*size(MyVals) &
       ,error_message('GatherAll: AllVals must be of size NumProc*size(MyVals)'))
    error_out = Epetra_MpiComm_GatherAll_Int(this%MpiComm_id,MyVals,AllVals,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%gather_int: failed.')
  end subroutine
 
  subroutine gather_long(this,MyVals,AllVals,err)
    class(Epetra_MpiComm)         ,intent(in)   :: this
    integer(c_long), dimension(:) ,intent(in)   :: MyVals
    integer(c_long), dimension(:) ,intent(inout):: AllVals
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    call assert(size(AllVals)==this%NumProc()*size(MyVals)&
       ,error_message('gather_long: AllVals must be of size NumProc*size(MyVals)'))
    error_out = Epetra_MpiComm_GatherAll_Long(this%MpiComm_id,MyVals,AllVals,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%gather_long: failed.')
  end subroutine

  subroutine sum_double(this,PartialSums,GlobalSums,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialSums 
    real(c_double), dimension(:) ,intent(inout) :: GlobalSums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_SumAll_Double(this%MpiComm_id,PartialSums,GlobalSums,size(PartialSums))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%sum_double: failed.')
  end subroutine
  
  subroutine sum_int(this,PartialSums,GlobalSums,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialSums 
    integer(c_int), dimension(:) ,intent(inout) :: GlobalSums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_SumAll_Int(this%MpiComm_id,PartialSums,GlobalSums,size(PartialSums))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%sum_int: failed.')
  end subroutine

  subroutine sum_long(this,PartialSums,GlobalSums,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialSums 
    integer(c_long), dimension(:),intent(inout) :: GlobalSums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_SumAll_Long(this%MpiComm_id,PartialSums,GlobalSums,size(PartialSums))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%sum_long: failed.')
  end subroutine

  subroutine max_double(this,PartialMaxs,GlobalMaxs,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialMaxs 
    real(c_double), dimension(:) ,intent(inout) :: GlobalMaxs
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_MaxAll_Double(this%MpiComm_id,PartialMaxs,GlobalMaxs,size(PartialMaxs))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%max_double: failed.')
  end subroutine

  subroutine max_int(this,PartialMaxs,GlobalMaxs,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialMaxs 
    integer(c_int), dimension(:) ,intent(inout) :: GlobalMaxs
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_MaxAll_Int(this%MpiComm_id,PartialMaxs,GlobalMaxs,size(PartialMaxs))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%max_int: failed.')
  end subroutine

  subroutine max_long(this,PartialMaxs,GlobalMaxs,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialMaxs 
    integer(c_long), dimension(:),intent(inout) :: GlobalMaxs
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_MaxAll_Long(this%MpiComm_id,PartialMaxs,GlobalMaxs,size(PartialMaxs))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%max_long: failed.')
  end subroutine
  
  subroutine min_double(this,PartialMins,GlobalMins,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialMins 
    real(c_double), dimension(:) ,intent(inout) :: GlobalMins
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_MinAll_Double(this%MpiComm_id,PartialMins,GlobalMins,size(PartialMins))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%min_double: failed.')
  end subroutine

  subroutine min_int(this,PartialMins,GlobalMins,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialMins 
    integer(c_int), dimension(:) ,intent(inout) :: GlobalMins
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_MinAll_Int(this%MpiComm_id,PartialMins,GlobalMins,size(PartialMins))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%min_int: failed.')
  end subroutine

  subroutine min_long(this,PartialMins,GlobalMins,err)
    class(Epetra_MpiComm)         ,intent(in)    :: this
    integer(c_long), dimension(:) ,intent(in)    :: PartialMins 
    integer(c_long), dimension(:) ,intent(inout) :: GlobalMins
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_MinAll_Long(this%MpiComm_id,PartialMins,GlobalMins,size(PartialMins))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%min_long: failed.')
  end subroutine
 
  subroutine ScanSum_double(this,MyVals,scan_sums,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: MyVals
    real(c_double), dimension(:) ,intent(inout) :: scan_sums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_ScanSum_Double(this%MpiComm_id,MyVals,scan_sums,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%ScanSum_Double: failed.')
  end subroutine

  subroutine ScanSum_int(this,MyVals,scan_sums,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: MyVals
    integer(c_int), dimension(:) ,intent(inout) :: scan_sums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_ScanSum_Int(this%MpiComm_id,MyVals,scan_sums,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%ScanSum_int: failed.')
  end subroutine

  subroutine ScanSum_long(this,MyVals,scan_sums,err)
    class(Epetra_MpiComm)        ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: MyVals
    integer(c_long), dimension(:),intent(inout) :: scan_sums
    type(error) ,optional, intent(inout) :: err
    integer(c_int)     :: error_out
    error_out = Epetra_MpiComm_ScanSum_Long(this%MpiComm_id,MyVals,scan_sums,size(MyVals))
    if (present(err)) err=error(error_out,'Epetra_MpiComm%ScanSum_long: failed.')
  end subroutine

  integer(c_int) function MyPID(this)
   class(Epetra_MpiComm)     , intent(in) :: this
   MyPID=Epetra_MpiComm_MyPID(this%MpiComm_id)
  end function

  integer(c_int) function NumProc(this)
   class(Epetra_MpiComm)     , intent(in) :: this
   NumProc=Epetra_MpiComm_NumProc(this%MpiComm_id)
  end function

  subroutine invalidate_EpetraMpiComm_ID(this)
    class(Epetra_MpiComm) ,intent(inout) :: this
    call this%invalidate_EpetraComm_ID
    this%MpiComm_id%table = FT_Invalid_ID
    this%MpiComm_id%index = FT_Invalid_Index 
    this%MpiComm_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraMpiComm(this)
    class(Epetra_MpiComm) ,intent(inout) :: this
    call this%ctrilinos_delete_EpetraComm()
    call Epetra_MpiComm_Destroy(this%MpiComm_id)
  end subroutine
#endif
end module 
