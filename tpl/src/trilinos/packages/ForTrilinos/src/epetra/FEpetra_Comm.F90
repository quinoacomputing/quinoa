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
module FEpetra_Comm
  use ForTrilinos_universal ,only : universal
  use ForTrilinos_enums !,only: FT_Epetra_Comm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_error
  use ForTrilinos_table_man
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: Epetra_Comm ! Expose type/methods

  type ,abstract ,extends(universal) :: Epetra_Comm
   private
   type(FT_Epetra_Comm_ID_t)         :: comm_id 
  contains
    !Developers only  -- to be called by developers from other ForTrilinos modules, not by end applications:
    procedure                                     :: invalidate_EpetraComm_ID
    procedure                                     :: ctrilinos_delete_EpetraComm
    procedure                                     :: get_EpetraComm_ID
    procedure                                     :: set_EpetraComm_ID
    procedure                 ,nopass             :: alias_EpetraComm_ID
    procedure ,non_overridable                    :: generalize_EpetraComm
    !Barrier Methods
    procedure(barrier_interface)          ,deferred  ::barrier
    !Broadcast Methods
    procedure(broadcast_double_interface) ,deferred  ::broadcast_double
    procedure(broadcast_int_interface)    ,deferred  ::broadcast_int
    procedure(broadcast_long_interface)   ,deferred  ::broadcast_long
    procedure(broadcast_char_interface)   ,deferred  ::broadcast_char
    generic :: broadcast=>broadcast_double,broadcast_int,broadcast_char
    !Gather Methods
    procedure(gather_double_interface),deferred  ::gather_double
    procedure(gather_int_interface)   ,deferred  ::gather_int
    procedure(gather_long_interface)  ,deferred  ::gather_long
    generic :: GatherAll=>gather_double,gather_int
    !Sum Methods
    procedure(sum_double_interface),deferred   ::sum_double
    procedure(sum_int_interface)   ,deferred   ::sum_int
    procedure(sum_long_interface)  ,deferred   ::sum_long
    generic :: SumAll=>sum_double,sum_int
    !Max/Min Methods
    procedure(max_double_interface) ,deferred   ::max_double
    procedure(max_int_interface)    ,deferred   ::max_int
    procedure(max_long_interface)   ,deferred   ::max_long
    generic :: MaxAll=>max_double,max_int
    procedure(min_double_interface) ,deferred   ::min_double
    procedure(min_int_interface)    ,deferred   ::min_int
    procedure(min_long_interface)               ,deferred   ::min_long
    generic :: MinAll=>min_double,min_int
    !Parallel Prefix Methods
    procedure(ScanSum_double_interface) ,deferred   ::ScanSum_double
    procedure(ScanSum_int_interface)    ,deferred   ::ScanSum_int
    procedure(ScanSum_long_interface)   ,deferred   ::ScanSum_long
    generic :: ScanSum=>ScanSum_double,ScanSum_int
    !Attribute Accessor Methods
    procedure(MyPID_interface)           ,deferred::MyPID
    procedure(NumProc_interface)         ,deferred::NumProc
  end type
  
  abstract interface

    subroutine barrier_interface(this) 
      import:: Epetra_Comm
      class(Epetra_Comm) ,intent(in)  :: this
    end subroutine
    subroutine broadcast_double_interface(this,MyVals,root,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      real(c_double) ,dimension(:) ,intent(inout) :: MyVals
      integer(c_int)               ,intent(in)    :: root
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine broadcast_int_interface(this,MyVals,root,err) 
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm,error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_int) ,dimension(:) ,intent(inout) :: MyVals
      integer(c_int)               ,intent(in)    :: root
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine broadcast_long_interface(this,MyVals,root,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm,error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_long),dimension(:) ,intent(inout) :: MyVals
      integer(c_int)               ,intent(in)    :: root
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine broadcast_char_interface(this,MyVals,root,err) 
      use iso_c_binding ,only: c_int,c_char
      import:: Epetra_Comm, error
      class(Epetra_Comm)                 ,intent(in)    :: this
      character(kind=c_char),dimension(:),intent(inout) :: MyVals
      integer(c_int)                     ,intent(in)    :: root
      type(error)   ,optional            ,intent(inout) :: err
    end subroutine
    subroutine gather_double_interface(this,MyVals,AllVals,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)         ,intent(in)    :: this
      real(c_double),dimension(:),intent(in)    :: MyVals
      real(c_double),dimension(:),intent(inout) :: AllVals
      type(error)   ,optional    ,intent(inout) :: err
    end subroutine
    subroutine gather_int_interface(this,MyVals,AllVals,err) 
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)         ,intent(in)   :: this
      integer(c_int),dimension(:),intent(in)   :: MyVals
      integer(c_int),dimension(:),intent(inout):: AllVals
      type(error)   ,optional    ,intent(inout) :: err
    end subroutine
    subroutine gather_long_interface(this,MyVals,AllVals,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_long),dimension(:) ,intent(in)    :: MyVals
      integer(c_long),dimension(:) ,intent(inout):: AllVals
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine sum_double_interface(this,PartialSums,GlobalSums,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      real(c_double),dimension(:)  ,intent(in)    :: PartialSums
      real(c_double),dimension(:)  ,intent(inout) :: GlobalSums
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine sum_int_interface(this,PartialSums,GlobalSums,err) 
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_int), dimension(:) ,intent(in)    :: PartialSums
      integer(c_int), dimension(:) ,intent(inout) :: GlobalSums
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine sum_long_interface(this,PartialSums,GlobalSums,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_long), dimension(:),intent(in)    :: PartialSums
      integer(c_long), dimension(:),intent(inout) :: GlobalSums
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine max_double_interface(this,PartialMaxs,GlobalMaxs,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      real(c_double), dimension(:) ,intent(in)    :: PartialMaxs
      real(c_double), dimension(:) ,intent(inout) :: GlobalMaxs
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine max_int_interface(this,PartialMaxs,GlobalMaxs,err) 
      use iso_c_binding ,only: c_int,c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_int), dimension(:) ,intent(in)    :: PartialMaxs
      integer(c_int), dimension(:) ,intent(inout) :: GlobalMaxs
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine max_long_interface(this,PartialMaxs,GlobalMaxs,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_long), dimension(:),intent(in)    :: PartialMaxs
      integer(c_long), dimension(:),intent(inout) :: GlobalMaxs
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine min_double_interface(this,PartialMins,GlobalMins,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      real(c_double), dimension(:) ,intent(in)    :: PartialMins
      real(c_double), dimension(:) ,intent(inout) :: GlobalMins
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine min_int_interface(this,PartialMins,GlobalMins,err) 
      use iso_c_binding ,only: c_int,c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_int), dimension(:) ,intent(in)    :: PartialMins
      integer(c_int), dimension(:) ,intent(inout) :: GlobalMins
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine min_long_interface(this,PartialMins,GlobalMins,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_long), dimension(:),intent(in)    :: PartialMins
      integer(c_long), dimension(:),intent(inout) :: GlobalMins
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine ScanSum_double_interface(this,MyVals,scan_sums,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      real(c_double), dimension(:) ,intent(in)    :: MyVals
      real(c_double), dimension(:) ,intent(inout) :: scan_sums 
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine ScanSum_int_interface(this,MyVals,scan_sums,err) 
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_int), dimension(:) ,intent(in)    :: MyVals
      integer(c_int), dimension(:) ,intent(inout) :: scan_sums 
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    subroutine ScanSum_long_interface(this,MyVals,scan_sums,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_long), dimension(:),intent(in)    :: MyVals
      integer(c_long), dimension(:),intent(inout) :: scan_sums 
      type(error)   ,optional      ,intent(inout) :: err
    end subroutine
    function MyPID_interface(this)
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm
      class(Epetra_Comm), intent(in) :: this
      integer(c_int) :: MyPID_interface
    end function
    function NumProc_interface(this)
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm
      class(Epetra_Comm), intent(in) :: this
      integer(c_int) :: NumProc_interface
    end function
  end interface

  contains
  ! Common type-bound procedures begin here: all ForTrilinos abstract base classes of the API have procedures analogous to these.
  ! The function set_EpetraComm_ID  must construct a struct id for the extended child class
  ! The functions should be called by developers only from other ForTrilinos modules, never by end applications:
  
 !----------------- Struct access ---------------------------------------------

  function get_EpetraComm_ID(this)
    class(Epetra_Comm) ,intent(in) :: this
    type(FT_Epetra_Comm_ID_t) :: get_EpetraComm_ID
    get_EpetraComm_ID = this%comm_id
  end function
  
  subroutine set_EpetraComm_ID(this,id)
    class(Epetra_Comm)        ,intent(inout) :: this
    type(FT_Epetra_Comm_ID_t) ,intent(in)    :: id 
    this%comm_id=id
  end subroutine 
  
  !----------------- Type casting ---------------------------------------------

  function alias_EpetraComm_ID(generic_id)
    use iso_c_binding, only : c_loc,c_int
    use ForTrilinos_table_man
    use ForTrilinos_enums
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status 
    type(error) :: ierr
    type(FT_Epetra_Comm_ID_t) :: alias_EpetraComm_ID
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Comm_ID),stat=status)
    ierr=error(status,'FEpetra_Comm:alias_EpetraComm_ID')
    call ierr%check_success()
    alias_EpetraComm_ID=degeneralize_EpetraComm(c_loc(alias_id))
  end function

  function generalize_EpetraComm(this)
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_Comm) ,intent(in) ,target :: this
   type(ForTrilinos_Universal_ID_t) :: generalize_EpetraComm
   generalize_EpetraComm = generalize_all( c_loc(this%comm_id) )
  end function
  
  function degeneralize_EpetraComm(generic_id) 
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)              ,value  :: generic_id
    type(FT_Epetra_Comm_ID_t),pointer:: local_ptr=>null()
    type(FT_Epetra_Comm_ID_t) :: degeneralize_EpetraComm
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraComm = local_ptr
  end function

  !__________ Garbage collection __________________________________________________

  subroutine invalidate_EpetraComm_ID(this)
    class(Epetra_Comm) ,intent(inout) :: this
    this%comm_id%table =FT_Invalid_ID
    this%comm_id%index =FT_Invalid_Index
    this%comm_id%is_const=FT_FALSE
  end subroutine
 
  subroutine ctrilinos_delete_EpetraComm(this)
    class(Epetra_Comm) ,intent(inout) :: this
    call Epetra_Comm_Destroy( this%comm_id )
  end subroutine

end module 
