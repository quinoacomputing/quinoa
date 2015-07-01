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


module FEpetra_CrsGraph
  use ForTrilinos_enums ,only: FT_Epetra_CrsGraph_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_enum_wrappers ,only: FT_Epetra_DataAccess_E_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_BlockMap  ,only: Epetra_BlockMap 
  use iso_c_binding     ,only: c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: Epetra_CrsGraph ! Expose type/constructors/methods

  type ,extends(universal)                    :: Epetra_CrsGraph 
    private
    type(FT_Epetra_CrsGraph_ID_t) :: CrsGraph_id 
  contains
     !Constructors
     procedure ,private :: duplicate_
     procedure ,private :: Create_ 
     procedure ,private :: Create_VarPerRow_
     !ForTrilinos Developers only -- not to be called in end-user application code.
     procedure ,private :: from_struct_
     generic :: Epetra_CrsGraph_ => from_struct_,duplicate_,Create_,Create_VarPerRow_
     procedure         :: invalidate_id => invalidate_EpetraCrsGraph_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraCrsGraph
     procedure         :: get_EpetraCrsGraph_ID 
     procedure ,nopass :: alias_EpetraCrsGraph_ID
     procedure         :: generalize 
  end type

   interface Epetra_CrsGraph ! constructors
     module procedure duplicate,from_struct, Create, Create_VarPerRow
   end interface

contains
  subroutine from_struct_(this,id)
    class(Epetra_CrsGraph) ,intent(out) :: this
    type(FT_Epetra_CrsGraph_ID_t) ,intent(in) :: id
    this%CrsGraph_id = id
    call this%register_self
  end subroutine

  type(Epetra_CrsGraph) function from_struct(id) result(new_Epetra_CrsGraph)
     type(FT_Epetra_CrsGraph_ID_t) ,intent(in) :: id
     call new_Epetra_CrsGraph%Epetra_CrsGraph_(id)
  end function

  subroutine Create_(this,CV,RowMap,NumIndicesPerRow,StaticProfile)
    use ForTrilinos_enums ,only: FT_boolean_t,FT_TRUE,FT_FALSE
    use iso_c_binding     ,only: c_int
    class(Epetra_CrsGraph) ,intent(out) :: this
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_BlockMap) ,intent(in) :: RowMap
    integer(c_int)         ,intent(in) :: NumIndicesPerRow
    logical,        optional                   :: StaticProfile        
    integer(FT_boolean_t)                      :: StaticProfile_in
    if (.not.present(StaticProfile)) then
      StaticProfile_in=FT_FALSE
    elseif (StaticProfile) then
      StaticProfile_in=FT_TRUE
    else
      StaticProfile_in=FT_FALSE
    endif
    call this%Epetra_CrsGraph_(Epetra_CrsGraph_Create(CV,RowMap%get_EpetraBlockMap_ID(),NumIndicesPerRow,StaticProfile_in))
  end subroutine

  function Create(CV,RowMap,NumIndicesPerRow,StaticProfile) result(new_Epetra_CrsGraph)
    use iso_c_binding     ,only: c_int
    type(Epetra_CrsGraph) :: new_Epetra_CrsGraph
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_BlockMap) ,intent(in) :: RowMap
    integer(c_int)         ,intent(in) :: NumIndicesPerRow
    logical,        optional                   :: StaticProfile        
    call new_Epetra_CrsGraph%Epetra_CrsGraph_(CV,RowMap,NumIndicesPerRow,StaticProfile)
  end function

  subroutine Create_VarPerRow_(this,CV,RowMap,NumIndicesPerRow,StaticProfile)
    use ForTrilinos_enums ,only: FT_boolean_t,FT_TRUE,FT_FALSE
    use iso_c_binding     ,only: c_int
    class(Epetra_CrsGraph) ,intent(out) :: this
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_BlockMap) ,intent(in) :: RowMap
    integer(c_int),dimension(:) ,intent(in) :: NumIndicesPerRow
    logical,        optional                   :: StaticProfile        
    integer(FT_boolean_t)                      :: StaticProfile_in
    type(FT_Epetra_CrsGraph_ID_t)              :: Create_VarPerRow_id
    if (.not.present(StaticProfile)) then
      StaticProfile_in=FT_FALSE
    elseif (StaticProfile) then
      StaticProfile_in=FT_TRUE
    else
      StaticProfile_in=FT_FALSE
    endif
    call this% &
    Epetra_CrsGraph_(Epetra_CrsGraph_Create_VarPerRow(CV,RowMap%get_EpetraBlockMap_ID(),NumIndicesPerRow,StaticProfile_in))
  end subroutine

  function Create_VarPerRow(CV,RowMap,NumIndicesPerRow,StaticProfile) result(new_Epetra_CrsGraph)
    type(Epetra_CrsGraph) :: new_Epetra_CrsGraph
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_BlockMap) ,intent(in) :: RowMap
    integer(c_int),dimension(:) ,intent(in) :: NumIndicesPerRow
    logical,        optional                   :: StaticProfile        
    if (present(StaticProfile)) then
      call new_Epetra_CrsGraph%Epetra_CrsGraph_(CV,RowMap,NumIndicesPerRow,StaticProfile)
    else
      call new_Epetra_CrsGraph%Epetra_CrsGraph_(CV,RowMap,NumIndicesPerRow)
    end if
  end function

  subroutine duplicate_(this,copy)
    class(Epetra_CrsGraph) ,intent(in) :: this
    type(Epetra_CrsGraph) ,intent(out) :: copy
    call copy%Epetra_CrsGraph_(Epetra_CrsGraph_Duplicate(this%CrsGraph_id))
  end subroutine

  type(Epetra_CrsGraph) function duplicate(original)
    type(Epetra_CrsGraph) ,intent(in) :: original
    call original%Epetra_CrsGraph_(duplicate)
  end function

  type(FT_Epetra_CrsGraph_ID_t) function get_EpetraCrsGraph_ID(this)
    class(Epetra_CrsGraph) ,intent(in) :: this 
    get_EpetraCrsGraph_ID=this%CrsGraph_id
  end function
  
  type(FT_Epetra_CrsGraph_ID_t) function alias_EpetraCrsGraph_ID(generic_id)
    use ForTrilinos_table_man,only: CT_Alias
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t,FT_Epetra_CrsGraph_ID
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_CrsGraph_ID),stat=status)
    ierr=error(status,'FEpetra_CrsGraph:alias_EpetraCrsGraph_ID')
    call ierr%check_success()
    alias_EpetraCrsGraph_ID=degeneralize_EpetraCrsGraph(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only : c_loc
   class(Epetra_CrsGraph) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%CrsGraph_id))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_CrsGraph) ,intent(in) ,target :: this
   ! generalize = Epetra_CrsGraph_Generalize ( this%CrsGraph_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_CrsGraph_ID_t) function degeneralize_EpetraCrsGraph(generic_id) bind(C)
  ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_CrsGraph_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                      ,value   :: generic_id
    type(FT_Epetra_CrsGraph_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraCrsGraph = local_ptr
  ! ____ Use for ForTrilinos function implementation ______

  ! ____ Use for CTrilinos function implementation ______
  !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
  !degeneralize_EpetraCrsGraph = Epetra_CrsGraph_Degeneralize(generic_id)
  ! ____ Use for CTrilinos function implementation ______
  end function

subroutine invalidate_EpetraCrsGraph_ID(this)
    class(Epetra_CrsGraph),intent(inout) :: this
    this%CrsGraph_id%table = FT_Invalid_ID
    this%CrsGraph_id%index = FT_Invalid_Index 
    this%CrsGraph_id%is_const = FT_FALSE
  end subroutine 

  subroutine ctrilinos_delete_EpetraCrsGraph(this)
    class(Epetra_CrsGraph),intent(inout) :: this
    call Epetra_CrsGraph_Destroy( this%CrsGraph_id ) 
  end subroutine

end module 

