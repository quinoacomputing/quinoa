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


module FEpetra_OffsetIndex
  use ForTrilinos_enums ,only: FT_Epetra_OffsetIndex_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_Import ,only: Epetra_Import
  use FEpetra_Export ,only: Epetra_Export
  use FEpetra_CrsGraph ,only: Epetra_CrsGraph
  use iso_c_binding     ,only: c_int,c_double,c_char
  use forepetra
  implicit none
  private                      ! Hide everything by default
  public :: Epetra_OffsetIndex ! Expose type/constructors/methods

  type ,extends(universal)                    :: Epetra_OffsetIndex 
    private
    type(FT_Epetra_OffsetIndex_ID_t) :: OffsetIndex_id 
  contains
     procedure ,private :: duplicate_
     procedure ,private :: create_FromExporter_
     procedure ,private :: create_FromImporter_
     procedure ,private :: from_struct_
     generic :: Epetra_OffsetIndex_ => duplicate_,from_struct_,create_FromExporter_,create_FromImporter_
     ! ForTrilinos developers only -- not for use in end-user application:
     procedure         :: invalidate_id => invalidate_EpetraOffsetIndex_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraOffsetIndex
     procedure         :: get_EpetraOffsetIndex_ID 
     procedure ,nopass :: alias_EpetraOffsetIndex_ID
     procedure         :: generalize 
  end type

   interface Epetra_OffsetIndex ! constructors
     module procedure duplicate,from_struct,create_FromExporter, create_FromImporter
   end interface

contains

  subroutine from_struct_(this,id)
    class(Epetra_OffsetIndex) ,intent(out) :: this
    type(FT_Epetra_OffsetIndex_ID_t) ,intent(in) :: id
    this%OffsetIndex_id = id
    call this%register_self
  end subroutine 

  function from_struct(id) result(new_Epetra_OffsetIndex)
    type(Epetra_OffsetIndex) :: new_Epetra_OffsetIndex
    type(FT_Epetra_OffsetIndex_ID_t) ,intent(in) :: id
    call new_Epetra_OffsetIndex%Epetra_OffsetIndex_(id)
  end function

  subroutine create_FromExporter_(this,SourceGraph,TargetGraph,exporter)
    class(Epetra_OffsetIndex) ,intent(out) :: this
    class(Epetra_CrsGraph),intent(in) :: SourceGraph,TargetGraph
    type(Epetra_Export), intent(in) :: exporter
    type(FT_Epetra_OffsetIndex_ID_t) :: create_FromExporter_id
    create_FromExporter_id = Epetra_OffsetIndex_Create_FromExporter(SourceGraph%get_EpetraCrsGraph_ID(),&
          TargetGraph%get_EpetraCrsGraph_ID(),exporter%get_EpetraExport_ID())
    call this%Epetra_OffsetIndex_(create_FromExporter_id)
  end subroutine

  function create_FromExporter(SourceGraph,TargetGraph,exporter) result(new_Epetra_OffsetIndex)
    type(Epetra_OffsetIndex) :: new_Epetra_OffsetIndex
    class(Epetra_CrsGraph),intent(in) :: SourceGraph,TargetGraph
    type(Epetra_Export), intent(in) :: exporter
    call new_Epetra_OffsetIndex%Epetra_OffsetIndex_(SourceGraph,TargetGraph,exporter)
  end function

  subroutine create_FromImporter_(this,SourceGraph,TargetGraph,importer)
    class(Epetra_OffsetIndex) ,intent(out) :: this
    class(Epetra_CrsGraph),intent(in) :: SourceGraph,TargetGraph
    type(Epetra_Import), intent(in) :: importer
    type(FT_Epetra_OffsetIndex_ID_t) :: create_FromImporter_id
    create_FromImporter_id = Epetra_OffsetIndex_Create_FromImporter(SourceGraph%get_EpetraCrsGraph_ID(),&
         TargetGraph%get_EpetraCrsGraph_ID(),importer%get_EpetraImport_ID())
    call this%Epetra_OffsetIndex_(create_FromImporter_id)
  end subroutine 

  function create_FromImporter(SourceGraph,TargetGraph,importer) result(new_Epetra_OffsetIndex)
    type(Epetra_OffsetIndex) :: new_Epetra_OffsetIndex
    class(Epetra_CrsGraph),intent(in) :: SourceGraph,TargetGraph
    type(Epetra_Import), intent(in) :: importer
    call new_Epetra_OffsetIndex%Epetra_OffsetIndex_(SourceGraph,TargetGraph,importer)
  end function

  subroutine duplicate_(this,copy)
    class(Epetra_OffsetIndex) ,intent(in) :: this
    type(Epetra_OffsetIndex) ,intent(out) :: copy
    call copy%Epetra_OffsetIndex_(Epetra_OffsetIndex_Duplicate(this%OffsetIndex_id))
  end subroutine

  type(Epetra_OffsetIndex) function duplicate(original)
    type(Epetra_OffsetIndex) ,intent(in) :: original
    call original%Epetra_OffsetIndex_(duplicate)
  end function

  type(FT_Epetra_OffsetIndex_ID_t) function get_EpetraOffsetIndex_ID(this)
    class(Epetra_OffsetIndex) ,intent(in) :: this 
    get_EpetraOffsetIndex_ID=this%OffsetIndex_id
  end function
  
  type(FT_Epetra_OffsetIndex_ID_t) function alias_EpetraOffsetIndex_ID(generic_id)
    use ForTrilinos_table_man,only: CT_Alias
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t,FT_Epetra_OffsetIndex_ID
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_OffsetIndex_ID),stat=status)
    ierr=error(status,'FEpetra_OffsetIndex:alias_EpetraOffsetIndex_ID')
    call ierr%check_success()
    alias_EpetraOffsetIndex_ID=degeneralize_EpetraOffsetIndex(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only : c_loc
   class(Epetra_OffsetIndex) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%OffsetIndex_id))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_OffsetIndex) ,intent(in) ,target :: this
   ! generalize = Epetra_OffsetIndex_Generalize ( this%OffsetIndex_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_OffsetIndex_ID_t) function degeneralize_EpetraOffsetIndex(generic_id) bind(C)
  ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_OffsetIndex_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                      ,value   :: generic_id
    type(FT_Epetra_OffsetIndex_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraOffsetIndex = local_ptr
  ! ____ Use for ForTrilinos function implementation ______

  ! ____ Use for CTrilinos function implementation ______
  !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
  !degeneralize_EpetraOffsetIndex = Epetra_OffsetIndex_Degeneralize(generic_id)
  ! ____ Use for CTrilinos function implementation ______
  end function

 subroutine invalidate_EpetraOffsetIndex_ID(this)
   class(Epetra_OffsetIndex),intent(inout) :: this  
   this%OffsetIndex_id%table = FT_Invalid_ID
   this%OffsetIndex_id%index = FT_Invalid_Index 
   this%OffsetIndex_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraOffsetIndex(this)
    class(Epetra_OffsetIndex),intent(inout) :: this
    call Epetra_OffsetIndex_Destroy( this%OffsetIndex_id ) 
  end subroutine

end module 
