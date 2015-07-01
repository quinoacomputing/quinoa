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
module FEpetra_DistObject
  use ForTrilinos_enums ,only : FT_Epetra_SrcDistObject_ID,FT_Epetra_DistObject_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
  use ForTrilinos_universal
  use ForTrilinos_table_man
  use ForTrilinos_error ,only : error
  use FEpetra_SrcDistObject ,only : Epetra_SrcDistObject
  use FEpetra_Export, only: Epetra_Export
  use FEpetra_Import, only: Epetra_Import
  use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
  use iso_c_binding     ,only : c_int,c_long,c_double,c_char
  use forepetra
  implicit none
  private                     ! Hide everything by default
  public :: Epetra_DistObject ! Expose type/constructors/methods

  type ,extends(Epetra_SrcDistObject)        :: Epetra_DistObject 
    private
    type(FT_Epetra_DistObject_ID_t) :: DistObject_id 
  contains
     !Constructor
     !Import/Export methods
     procedure, private:: DistObject_Export
     procedure, private:: DistObject_Export_UsingImporter
     generic :: export => DistObject_Export_UsingImporter, DistObject_Export
     procedure, private:: DistObject_Import
     procedure, private:: DistObject_Import_UsingExporter
     generic :: import => DistObject_Import_UsingExporter, DistObject_Import
     !Developers only
     procedure ,private :: from_struct__
     generic :: Epetra_DistObject_ => from_struct__
     procedure         :: invalidate_id => invalidate_EpetraDistObject_ID 
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraDistObject
     procedure         :: get_EpetraDistObject_ID 
     procedure ,nopass :: alias_EpetraDistObject_ID
     procedure         :: generalize 
  end type

   interface Epetra_DistObject ! constructors
     module procedure from_struct
   end interface

contains

  subroutine from_struct__(this,id)
    class(Epetra_DistObject) ,intent(out) :: this
    type(FT_Epetra_DistObject_ID_t) ,intent(in) :: id
    this%DistObject_id = id
    call this%set_EpetraSrcDistObject_ID(this%alias_EpetraSrcDistObject_ID(this%generalize()))
    call this%register_self
  end subroutine

  function from_struct(id) result(new_Epetra_DistObject)
    type(Epetra_DistObject) :: new_Epetra_DistObject
    type(FT_Epetra_DistObject_ID_t) ,intent(in) :: id
    call new_Epetra_DistObject%Epetra_DistObject_(id)
  end function

  type(FT_Epetra_DistObject_ID_t) function get_EpetraDistObject_ID(this)
   class(Epetra_DistObject) ,intent(in) :: this 
   get_EpetraDistObject_ID=this%DistObject_id
  end function
  
  type(FT_Epetra_DistObject_ID_t) function alias_EpetraDistObject_ID(generic_id)
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: FT_Epetra_DistObject_ID,ForTrilinos_Universal_ID_t
    use ForTrilinos_table_man,only: CT_Alias
    type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(Fortrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_DistObject_ID),stat=status)
    ierr=error(status,'FEpetra_DistObject:alias_EpetraDistObject_ID')
    call ierr%check_success()
    alias_EpetraDistObject_ID=degeneralize_EpetraDistObject(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_DistObject) ,intent(in) ,target :: this
   generalize = generalize_all( c_loc(this%DistObject_id) )
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_DistObject) ,intent(in) ,target :: this
   ! generalize = Epetra_DistObject_Generalize ( this%DistObject_id ) 
   ! ____ Use for CTrilinos function implementation ______
  end function
 
 type(FT_Epetra_DistObject_ID_t) function degeneralize_EpetraDistObject(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_DistObject_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                     ,value   :: generic_id
    type(FT_Epetra_DistObject_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraDistObject = local_ptr
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
   ! degeneralize_EpetraDistObject = Epetra_DistObject_Degeneralize( generic_id )
   ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine DistObject_Export(this,A,exporter,CombineMode,indexor,err)
   class(Epetra_DistObject), intent(in) :: this 
   class(Epetra_SrcDistObject), intent(in) :: A
   type(Epetra_Export),intent(in) :: exporter
   integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
   type(Epetra_OffsetIndex), intent(in) :: indexor
   type(error),optional,intent(out) :: err
   integer(c_int)     :: error_out
   error_out=Epetra_DistObject_Export(this%DistObject_id,A%get_EpetraSrcDistObject_ID(),&
                exporter%get_EpetraExport_ID(),CombineMode,indexor%get_EpetraOffsetIndex_ID())
   if (present(err)) err=error(error_out,'Epetra_DistObject%Export: failed.')
  end subroutine

  subroutine DistObject_Export_UsingImporter(this,A,importer,CombineMode,indexor,err)
   class(Epetra_DistObject), intent(in) :: this
   class(Epetra_SrcDistObject), intent(in) :: A
   type(Epetra_Import),intent(in) :: importer
   integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
   type(Epetra_OffsetIndex), intent(in) :: indexor
   type(error),optional,intent(out) :: err
   integer(c_int)     :: error_out
   error_out=Epetra_DistObject_Export_UsingImporter(this%DistObject_id,A%get_EpetraSrcDistObject_ID(),&
          importer%get_EpetraImport_ID(),CombineMode,indexor%get_EpetraOffsetIndex_ID())
   if (present(err)) err=error(error_out,'Epetra_DistObject%Export_UsingImporter: failed.')
  end subroutine

  subroutine DistObject_Import(this,A,importer,CombineMode,indexor,err)
   class(Epetra_DistObject), intent(in) :: this
   class(Epetra_SrcDistObject), intent(in) :: A
   type(Epetra_Import),intent(in) :: importer
   integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
   type(Epetra_OffsetIndex), intent(in) :: indexor
   type(error),optional,intent(out) :: err
   integer(c_int)     :: error_out
   error_out=Epetra_DistObject_Import(this%DistObject_id,A%get_EpetraSrcDistObject_ID(),&
       importer%get_EpetraImport_ID(),CombineMode,indexor%get_EpetraOffsetIndex_ID())
   if (present(err)) err=error(error_out,'Epetra_DistObject%Import: failed.')
  end subroutine

  subroutine DistObject_Import_UsingExporter(this,A,exporter,CombineMode,indexor,err)
   class(Epetra_DistObject), intent(in) :: this
   class(Epetra_SrcDistObject), intent(in) :: A
   type(Epetra_Export),intent(in) :: exporter
   integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
   type(Epetra_OffsetIndex), intent(in) :: indexor
   type(error),optional,intent(out) :: err
   integer(c_int)     :: error_out
   error_out=Epetra_DistObject_Import_UsingExporter(this%DistObject_id,A%get_EpetraSrcDistObject_ID(),&
         exporter%get_EpetraExport_ID(),CombineMode,indexor%get_EpetraOffsetIndex_ID())
   if (present(err)) err=error(error_out,'Epetra_DistObject%Import: failed.')
  end subroutine

  subroutine invalidate_EpetraDistObject_ID(this)
    class(Epetra_DistObject),intent(inout) :: this
    call this%invalidate_EpetraSrcDistObject_ID
    this%DistObject_id%table=FT_Invalid_ID
    this%DistObject_id%index=FT_Invalid_Index
    this%DistObject_id%is_const=FT_FALSE
  end subroutine
  
  subroutine ctrilinos_delete_EpetraDistObject(this)
    class(Epetra_DistObject) ,intent(inout) :: this
    call this%ctrilinos_delete_EpetraSrcDistObject()
    call Epetra_DistObject_Destroy(this%DistObject_id)
  end subroutine

end module 
