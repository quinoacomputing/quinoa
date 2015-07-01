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


module FEpetra_Import
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_Import_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal ,only:universal
  use ForTrilinos_error
  use FEpetra_Comm  ,only: Epetra_Comm
  use FEpetra_BlockMap ,only: Epetra_BlockMap
  use iso_c_binding ,only: c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: Epetra_Import ! Expose type/constructors/methods

  type ,extends(universal) :: Epetra_Import 
    private
    type(FT_Epetra_Import_ID_t) :: Import_id 
  contains
     ! Constructors
     procedure ,private :: create_
     procedure ,private :: duplicate_
     procedure ,private :: from_struct_
     generic :: Epetra_Import_ => create_,duplicate_,from_struct_
     ! Public member functions
     procedure        :: NumSameIDs
     procedure        :: NumPermuteIDs
     procedure        :: PermuteFromLIDs
     procedure        :: PermuteToLIDs
     procedure        :: NumRemoteIDs
     procedure        :: NumExportIDs
     procedure        :: NumSend
     procedure        :: NumRecv
     procedure        :: SourceMap
     procedure        :: TargetMap
     !Developers only
     procedure         :: invalidate_id => invalidate_EpetraImport_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraImport
     procedure         :: get_EpetraImport_ID 
     procedure ,nopass :: alias_EpetraImport_ID
     procedure         :: generalize 
  end type

   interface Epetra_Import ! constructors
     module procedure create,duplicate,from_struct
   end interface
 
contains

  subroutine from_struct_(this,id)
    class(Epetra_Import) ,intent(out) :: this
    type(FT_Epetra_Import_ID_t) ,intent(in) :: id
    this%Import_id = id
    call this%register_self 
  end subroutine
 
  function from_struct(id) result(new_Epetra_Import)
    type(Epetra_Import) :: new_Epetra_Import
    type(FT_Epetra_Import_ID_t) ,intent(in) :: id
    call new_Epetra_Import%Epetra_Import_(id)
  end function
 
  ! Original C++ prototype:
  ! Epetra_Import( const Epetra_BlockMap & TargetMap, const Epetra_BlockMap & SourceMap );
  ! CTrilinos prototype:
  ! CT_Epetra_Import_ID_t Epetra_Import_Create ( CT_Epetra_BlockMap_ID_t TargetMapID,
  ! CT_Epetra_BlockMap_ID_t SourceMapID );

  subroutine create_(this,TargetMap,SourceMap)
    class(Epetra_Import) ,intent(out) :: this
    class(Epetra_BlockMap), intent(in) :: TargetMap,SourceMap
    call this%Epetra_Import_(Epetra_Import_Create(TargetMap%get_EpetraBlockMap_ID(),SourceMap%get_EpetraBlockMap_ID()))
  end subroutine

  function create(TargetMap,SourceMap) result(new_Epetra_Import)
    type(Epetra_Import) :: new_Epetra_Import
    class(Epetra_BlockMap), intent(in) :: TargetMap,SourceMap
    call new_Epetra_Import%Epetra_Import_(TargetMap,SourceMap)
  end function

  ! Original C++ prototype:
  ! Epetra_Import(const Epetra_Import& Importer);
  ! CTrilinos prototype:
  ! CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( CT_Epetra_Import_ID_t ImporterID );

  subroutine duplicate_(this,copy)
    class(Epetra_Import) ,intent(in) :: this
    type(Epetra_Import) ,intent(out) :: copy
    call copy%Epetra_Import_(Epetra_Import_Duplicate(this%import_id))
  end subroutine

  type(Epetra_Import) function duplicate(original)
    type(Epetra_Import) ,intent(in) :: original
    call original%Epetra_Import_(duplicate)
  end function


  type(FT_Epetra_Import_ID_t) function get_EpetraImport_ID(this)
    class(Epetra_Import) ,intent(in) :: this 
    get_EpetraImport_ID=this%Import_id
  end function
  
  type(FT_Epetra_Import_ID_t) function alias_EpetraImport_ID(generic_id)
    use ForTrilinos_table_man
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Epetra_Import_ID
    use iso_c_binding     ,only: c_loc,c_int
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Import_ID),stat=status)
    ierr=error(status,'FEpetra_Import:alias_EpetraImport_ID')
    call ierr%check_success()
    alias_EpetraImport_ID=degeneralize_EpetraImport(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(Epetra_Import) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%Import_ID))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(Epetra_Import) ,intent(in) ,target :: this
   !generalize = Epetra_Import_Generalize ( this%Import_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_Import_ID_t) function degeneralize_EpetraImport(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_Import_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                   ,value   :: generic_id
    type(FT_Epetra_Import_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraImport = local_ptr
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
   !degeneralize_EpetraImport = Epetra_Import_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  integer(c_int) function NumSameIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumSameIDs=Epetra_Import_NumSameIDs(this%Import_id)
  end function

  integer(c_int) function NumPermuteIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumPermuteIDs=Epetra_Import_NumPermuteIDs(this%Import_id)
  end function

 function PermuteFromLIDs(this)
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer,c_int
    class(Epetra_Import), intent(in) :: this
    integer(c_int),dimension(:),allocatable :: PermuteFromLIDs
    type(c_ptr)   :: PermuteFromLIDs_external_ptr 
    integer(c_int),pointer :: PermuteFromLIDs_local_ptr=>null()
    integer(c_int) :: status
    type(error) :: ierr
    allocate(PermuteFromLIDs(this%NumPermuteIDs()),stat=status)
    ierr=error(status,'FEpetra_Import:alias_EpetraImport_ID')
    call ierr%check_success()
    PermuteFromLIDs_external_ptr=Epetra_Import_PermuteFromLIDs(this%Import_id)
    call c_f_pointer (PermuteFromLIDs_external_ptr, PermuteFromLIDs_local_ptr)
    PermuteFromLIDs=PermuteFromLIDs_local_ptr
  end function

 function PermuteToLIDs(this)
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer,c_int
    class(Epetra_Import), intent(in) :: this
    integer(c_int),dimension(:),allocatable :: PermuteToLIDs
    type(c_ptr)   :: PermuteToLIDs_external_ptr 
    integer(c_int),pointer :: PermuteToLIDs_local_ptr=>null()
    integer(c_int) :: status
    type(error) :: ierr
    allocate(PermuteToLIDs(this%NumPermuteIDs()),stat=status)
    ierr=error(status,'FEpetra_Import:alias_EpetraImport_ID')
    call ierr%check_success()
    PermuteToLIDs_external_ptr=Epetra_Import_PermuteToLIDs(this%Import_id)
    call c_f_pointer (PermuteToLIDs_external_ptr, PermuteToLIDs_local_ptr)
    PermuteToLIDs=PermuteToLIDs_local_ptr
  end function

  integer(c_int) function NumRemoteIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumRemoteIDs=Epetra_Import_NumRemoteIDs(this%Import_id)
  end function

  integer(c_int) function NumExportIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumExportIDs=Epetra_Import_NumExportIDs(this%Import_id)
  end function

  integer(c_int) function NumSend(this)
    class(Epetra_Import), intent(in) :: this
    NumSend=Epetra_Import_NumSend(this%Import_id)
  end function

  integer(c_int) function NumRecv(this)
    class(Epetra_Import), intent(in) :: this
    NumRecv=Epetra_Import_NumRecv(this%Import_id)
  end function

  type(Epetra_BlockMap) function SourceMap(this)
   class(Epetra_Import), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: SourceMap_id
   SourceMap_id=Epetra_Import_SourceMap(this%Import_id)
   SourceMap=Epetra_BlockMap(SourceMap_id)
  end function

  type(Epetra_BlockMap) function TargetMap(this)
   class(Epetra_Import), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: TargetMap_id
   TargetMap_id=Epetra_Import_TargetMap(this%Import_id)
   TargetMap=Epetra_BlockMap(TargetMap_id)
  end function

  subroutine invalidate_EpetraImport_ID(this)
    class(Epetra_Import),intent(inout) :: this
    this%Import_id%table = FT_Invalid_ID
    this%Import_id%index = FT_Invalid_Index 
    this%Import_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraImport(this)
    class(Epetra_Import),intent(inout) :: this
    call Epetra_Import_Destroy( this%Import_id ) 
  end subroutine
end module 
