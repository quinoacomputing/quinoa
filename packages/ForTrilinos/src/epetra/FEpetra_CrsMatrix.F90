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

module FEpetra_CrsMatrix
  use ForTrilinos_enums !,only : FT_Epetra_RowMatrix_ID,FT_Epetra_CrsMatrix_ID_t,ForTrilinos_Universal_ID_t,
                        !        FT_Epetra_Map_ID_t,FT_boolean_t,FT_FALSE,FT_TRUE
  use ForTrilinos_enum_wrappers
  use ForTrilinos_table_man
  use ForTrilinos_error
  use ForTrilinos_assertion_utility
  use FEpetra_RowMatrix ,only : Epetra_RowMatrix
  use FEpetra_Map     !  ,only : Epetra_Map
  use iso_c_binding     ,only : c_int,c_double
  use forepetra
  implicit none
  private                         ! Hide everything by default
  public :: Epetra_CrsMatrix ! Expose type/constructors/methods

  type ,extends(Epetra_RowMatrix)  :: Epetra_CrsMatrix 
    private
    type(FT_Epetra_CrsMatrix_ID_t) :: CrsMatrix_id 
  contains
    !Constructors
    procedure ,private :: Create_VarPerRow__
    procedure ,private :: Create__
    procedure ,private :: Create_VarPerRow_WithColMap__
    procedure ,private :: Create_WithColMap__
    procedure ,private :: duplicate__
    procedure ,private :: from_struct__
    generic :: Epetra_CrsMatrix_=> from_struct__,Create_VarPerRow__,&
         & Create__,Create_VarPerRow_WithColMap__,&
         & Create_WithColMap__,duplicate__
    !Insertion/Replace/SumInto methods
    procedure         :: PutScalar
    procedure         :: Scale
    procedure         :: InsertGlobalValues
    procedure         :: ReplaceGlobalValues
    !Transformation methods
    procedure,private :: FillComplete_Op
    procedure,private :: FillComplete_Map
    generic           :: FillComplete=>FillComplete_Op, FillComplete_Map 
    !Matrix data extraction routines
    procedure         :: ExtractGlobalRowCopy
    procedure         :: NumMyRowEntries
    procedure         :: NumMyEntries
    procedure         :: MaxNumEntries
    !Computational Methods
    procedure         :: Multiply_Vector
    procedure         :: Multiply => Multiply_MultiVector
    !Attribute Accessor Methods
    procedure         :: RowMatrixRowMap 
    procedure         :: RowMap
    procedure         :: NumGlobalEntries
    !Local/Global ID method
    procedure         :: MyGlobalRow
    !Developers only
    procedure         :: invalidate_id => invalidate_EpetraCrsMatrix_ID
    procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraCrsMatrix
    procedure         :: get_EpetraCrsMatrix_ID 
    procedure ,nopass :: alias_EpetraCrsMatrix_ID
    procedure         :: generalize 
  end type

   interface Epetra_CrsMatrix ! constructors
     module procedure Create_VarPerRow,Create,Create_VarPerRow_WithColMap,Create_WithColMap,duplicate,from_struct
   end interface

contains

  subroutine from_struct__(this,id)
    class(Epetra_CrsMatrix) ,intent(out) :: this
    type(FT_Epetra_CrsMatrix_ID_t) ,intent(in) :: id
    this%CrsMatrix_id = id
    call this%set_EpetraRowMatrix_ID(this%alias_EpetraRowMatrix_ID(this%generalize()))
    call this%register_self
  end subroutine
 
  function from_struct(id) result(new_Epetra_CrsMatrix)
    type(Epetra_CrsMatrix) :: new_Epetra_CrsMatrix
    type(FT_Epetra_CrsMatrix_ID_t) ,intent(in) :: id
    call new_Epetra_CrsMatrix%Epetra_CrsMatrix_(id)
  end function
 
  ! Original C++ prototype:
  ! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const int* NumEntriesPerRow,
  ! bool StaticProfile = false);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( CT_Epetra_DataAccess_E_t CV,
  !CT_Epetra_Map_ID_t RowMapID, const int * NumEntriesPerRow, boolean StaticProfile);

  subroutine Create_VarPerRow__(this,CV,Row_Map,NumEntriesPerRow,StaticProfile)
    use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_CrsMatrix) ,intent(out) :: this
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_Map),              intent(in) :: Row_Map
    integer(c_int), dimension(:),   intent(in) :: NumEntriesPerRow 
    logical,        optional                   :: StaticProfile                  
    integer(FT_boolean_t)                      :: StaticProfile_in
    if (.not.present(StaticProfile)) then
      StaticProfile_in=FT_FALSE
    elseif (StaticProfile) then
      StaticProfile_in=FT_TRUE
    else
      StaticProfile_in=FT_FALSE
    endif
    call this%Epetra_CrsMatrix_(Epetra_CrsMatrix_Create_VarPerRow(CV,Row_Map%get_EpetraMap_ID(),NumEntriesPerRow,StaticProfile_in))
  end subroutine

  function Create_VarPerRow(CV,Row_Map,NumEntriesPerRow,StaticProfile) result(new_Epetra_CrsMatrix)
    type(Epetra_CrsMatrix) :: new_Epetra_CrsMatrix
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_Map),              intent(in) :: Row_Map
    integer(c_int), dimension(:),   intent(in) :: NumEntriesPerRow 
    logical,        optional                   :: StaticProfile                  
    if (present(StaticProfile)) then
      call new_Epetra_CrsMatrix%Epetra_CrsMatrix_(CV,Row_Map,NumEntriesPerRow,StaticProfile)
    else
      call new_Epetra_CrsMatrix%Epetra_CrsMatrix_(CV,Row_Map,NumEntriesPerRow)
    endif
  end function

  ! Original C++ prototype:
  ! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int NumEntriesPerRow, bool StaticProfile = false);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  ! int NumEntriesPerRow, boolean StaticProfile );

  subroutine Create__(this,CV,Row_Map,NumEntriesPerRow,StaticProfile)
    use ForTrilinos_enums,         only: FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_CrsMatrix) ,intent(out) :: this
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_Map),              intent(in) :: Row_Map
    integer(c_int),                 intent(in) :: NumEntriesPerRow
    logical,        optional                   :: StaticProfile
    integer(FT_boolean_t)                      :: StaticProfile_in
    if (.not.present(StaticProfile)) then
      StaticProfile_in=FT_FALSE
    elseif (StaticProfile) then
      StaticProfile_in=FT_TRUE
    else
      StaticProfile_in=FT_FALSE
    endif
    call this%Epetra_CrsMatrix_(Epetra_CrsMatrix_Create(CV,Row_Map%get_EpetraMap_ID(),NumEntriesPerRow,StaticProfile_in))
  end subroutine

  function Create(CV,Row_Map,NumEntriesPerRow,StaticProfile) result(new_Epetra_CrsMatrix)
    type(Epetra_CrsMatrix) :: new_Epetra_CrsMatrix
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_Map),              intent(in) :: Row_Map
    integer(c_int),                 intent(in) :: NumEntriesPerRow
    logical,        optional                   :: StaticProfile
    if (present(StaticProfile)) then
      call new_Epetra_CrsMatrix%Epetra_CrsMatrix_(CV,Row_Map,NumEntriesPerRow,StaticProfile)
    else
      call new_Epetra_CrsMatrix%Epetra_CrsMatrix_(CV,Row_Map,NumEntriesPerRow)
    endif
  end function

  ! Original C++ prototype:
  ! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, const int* NumEntriesPerRow, 
  ! bool StaticProfile = false);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( CT_Epetra_DataAccess_E_t CV, 
  !            CT_Epetra_Map_ID_t RowMapID, CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, boolean StaticProfile );

  subroutine Create_VarPerRow_WithColMap__(this,CV,Row_Map,Col_Map,NumEntriesPerRow,StaticProfile)
    use ForTrilinos_enums,         only: FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_CrsMatrix) ,intent(out) :: this
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_Map),              intent(in) :: Row_Map,Col_Map
    integer(c_int), dimension(:) ,  intent(in) :: NumEntriesPerRow
    logical,        optional                   :: StaticProfile
    integer(FT_boolean_t)                      :: StaticProfile_in   
    type(FT_Epetra_CrsMatrix_ID_t):: Create_id
    if (.not.present(StaticProfile)) then
      StaticProfile_in=FT_FALSE
    elseif (StaticProfile) then
      StaticProfile_in=FT_TRUE
    else
      StaticProfile_in=FT_FALSE
    endif
    Create_id = Epetra_CrsMatrix_Create_VarPerRow_WithColMap(CV,Row_Map%get_EpetraMap_ID(),Col_Map%get_EpetraMap_ID(),&
                                                             NumEntriesPerRow,StaticProfile_in)
    call this%Epetra_CrsMatrix_(Create_id)
  end subroutine

  function Create_VarPerRow_WithColMap(CV,Row_Map,Col_Map,NumEntriesPerRow,StaticProfile) result(new_Epetra_CrsMatrix)
    type(Epetra_CrsMatrix) :: new_Epetra_CrsMatrix
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_Map),              intent(in) :: Row_Map,Col_Map
    integer(c_int), dimension(:) ,  intent(in) :: NumEntriesPerRow
    logical,        optional                   :: StaticProfile
    if (present(StaticProfile)) then
      call new_Epetra_CrsMatrix%Epetra_CrsMatrix_(CV,Row_Map,Col_Map,NumEntriesPerRow,StaticProfile)
    else
      call new_Epetra_CrsMatrix%Epetra_CrsMatrix_(CV,Row_Map,Col_Map,NumEntriesPerRow)
    endif
  end function

  !Original C++ prototype:
  ! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, int NumEntriesPerRow, 
  !                  bool StaticProfile = false);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  !                  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, boolean StaticProfile );

  subroutine Create_WithColMap__(this,CV,Row_Map,Col_Map,NumEntriesPerRow,StaticProfile)
    use ForTrilinos_enums,         only: FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_CrsMatrix) ,intent(out) :: this
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_Map),              intent(in) :: Row_Map,Col_Map
    integer(c_int),                 intent(in) :: NumEntriesPerRow
    logical,        optional                   :: StaticProfile
    integer(FT_boolean_t)                      :: StaticProfile_in
    type(FT_Epetra_CrsMatrix_ID_t):: Create_id
    if (.not.present(StaticProfile)) then
      StaticProfile_in=FT_FALSE
    elseif (StaticProfile) then
      StaticProfile_in=FT_TRUE
    else
      StaticProfile_in=FT_FALSE
    endif
    Create_id = Epetra_CrsMatrix_Create_WithColMap(CV,Row_Map%get_EpetraMap_ID(),Col_Map%get_EpetraMap_ID(),&
                                                   NumEntriesPerRow,StaticProfile_in)
    call this%Epetra_CrsMatrix_(Create_id)
  end subroutine

  function Create_WithColMap(CV,Row_Map,Col_Map,NumEntriesPerRow,StaticProfile) result(new_Epetra_CrsMatrix)
    type(Epetra_CrsMatrix) :: new_Epetra_CrsMatrix
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_Map),              intent(in) :: Row_Map,Col_Map
    integer(c_int),                 intent(in) :: NumEntriesPerRow
    logical,        optional                   :: StaticProfile
    if (present(StaticProfile)) then
      call new_Epetra_CrsMatrix%Epetra_CrsMatrix_(CV,Row_Map,Col_Map,NumEntriesPerRow,StaticProfile)
    else
      call new_Epetra_CrsMatrix%Epetra_CrsMatrix_(CV,Row_Map,Col_Map,NumEntriesPerRow)
    endif
  end function

  ! Original C++ prototype:
  ! Epetra_CrsMatrix(const Epetra_CrsMatrix& RowMatrix);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( CT_Epetra_BasicRowMatrix_ID_t RowMatrixID );

  type(FT_Epetra_CrsMatrix_ID_t) function get_EpetraCrsMatrix_ID(this)
   class(Epetra_CrsMatrix) ,intent(in) :: this 
   get_EpetraCrsMatrix_ID=this%CrsMatrix_id
  end function

  subroutine duplicate__(this,copy)
    class(Epetra_CrsMatrix) ,intent(in) :: this
    type(Epetra_CrsMatrix) ,intent(out) :: copy
    call copy%Epetra_CrsMatrix_(Epetra_CrsMatrix_Duplicate(this%CrsMatrix_id))
  end subroutine

  type(Epetra_CrsMatrix) function duplicate(original)
    type(Epetra_CrsMatrix) ,intent(in) :: original
    call original%Epetra_CrsMatrix_(duplicate)
  end function
  
  type(FT_Epetra_CrsMatrix_ID_t) function alias_EpetraCrsMatrix_ID(generic_id)
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: FT_Epetra_CrsMatrix_ID,ForTrilinos_Universal_ID_t
    use ForTrilinos_table_man,only: CT_Alias
    type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(Fortrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_CrsMatrix_ID),stat=status)
    ierr=error(status,'FEpetra_CrsMatrix:alias_EpetraCrsMatrix_ID')
    call ierr%check_success()
    alias_EpetraCrsMatrix_ID=degeneralize_EpetraCrsMatrix(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_CrsMatrix) ,intent(in) ,target :: this
   generalize = generalize_all( c_loc(this%CrsMatrix_id) )
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_CrsMatrix) ,intent(in) ,target :: this
   ! generalize = Epetra_CrsMatrix_Generalize ( this%CrsMatrix_id ) 
   ! ____ Use for CTrilinos function implementation ______
  end function
 
 type(FT_Epetra_CrsMatrix_ID_t) function degeneralize_EpetraCrsMatrix(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_CrsMatrix_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                     ,value   :: generic_id
    type(FT_Epetra_CrsMatrix_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraCrsMatrix = local_ptr
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
   ! degeneralize_EpetraCrsMatrix = Epetra_CrsMatrix_Degeneralize( generic_id )
   ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine PutScalar(this,scalar,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   real(c_double)           ,intent(in):: scalar
   type(error), optional, intent(out) :: err
   integer(c_int)                          :: error_out
   error_out=Epetra_CrsMatrix_PutScalar(this%CrsMatrix_id,scalar)
   if (present(err)) err=error(error_out,'Epetra_CrsMatrix%PutScalar: failed.')
  end subroutine
 
  subroutine Scale(this,ScalarConstant,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   real(c_double)           ,intent(in):: ScalarConstant
   type(error), optional  ,intent(out) :: err
   integer(c_int)                      :: error_out
   error_out=Epetra_CrsMatrix_Scale(this%CrsMatrix_id,ScalarConstant)
   if (present(err)) err=error(error_out,'Epetra_CrsMatrix%Scale: failed.')
  end subroutine

  subroutine InsertGlobalValues(this,GlobalRow,values,indices,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: GlobalRow
   real(c_double),dimension(:),intent(in):: values 
   integer(c_int),dimension(:),intent(in):: indices 
   type(error), optional, intent(out) :: err
   integer(c_int)                          :: error_out
   call assert(size(values)==size(indices) &
        ,error_message('InsertGlobalValues: values and indices should have the same size'))
   error_out=Epetra_CrsMatrix_InsertGlobalValues(this%CrsMatrix_id,GlobalRow,size(values),values,indices)
   if (present(err)) err=error(error_out,'Epetra_CrsMatrix%InsertGlobalValues')
  end subroutine

  subroutine ReplaceGlobalValues(this,GlobalRow,values,indices,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: GlobalRow
   real(c_double), dimension(:)        :: values 
   integer(c_int),    dimension(:)        :: indices 
   type(error), optional, intent(out) :: err
   integer(c_int)                          :: error_out
   call assert(size(values)==size(indices) &
        ,error_message('ReplaceGlobalValues: values and indices should have the same size'))
   error_out=Epetra_CrsMatrix_ReplaceGlobalValues(this%CrsMatrix_id,GlobalRow,size(values),values,indices)
   if (present(err)) err=error(error_out,'Epetra_CrsMatrix%ReplaceGlobalValues: failed.')
  end subroutine
  
  subroutine FillComplete_Op(this,OptimizeDataStorage,err)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   class(Epetra_CrsMatrix), intent(in) :: this
   logical,optional,        intent(in) :: OptimizeDataStorage
   type(error),optional,intent(out)    :: err
   integer(c_int)                      :: error_out
   integer(FT_boolean_t)               :: OptimizeDataStorage_in
   if (.not.present(OptimizeDataStorage)) then
     OptimizeDataStorage_in=FT_TRUE
   elseif (OptimizeDataStorage) then
     OptimizeDataStorage_in=FT_TRUE
   else
     OptimizeDataStorage_in=FT_FALSE
   endif
   error_out=Epetra_CrsMatrix_FillComplete(this%CrsMatrix_id,OptimizeDataStorage_in)
   if (present(err)) err=error(error_out,'Epetra_CrsMatrix%FillComplete_On: failed.')
  end subroutine

  subroutine FillComplete_Map(this,DomainMap,RangeMap,OptimizeDataStorage,err)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   class(Epetra_CrsMatrix), intent(in) :: this
   class(Epetra_Map) , intent(in) :: DomainMap
   class(Epetra_Map) , intent(in) :: RangeMap
   logical,optional,        intent(in) :: OptimizeDataStorage
   type(error),optional,intent(out)    :: err
   integer(c_int)                      :: error_out 
   integer(FT_boolean_t)               :: OptimizeDataStorage_in
   if (.not.present(OptimizeDataStorage)) then
     OptimizeDataStorage_in=FT_TRUE
   elseif (OptimizeDataStorage) then
     OptimizeDataStorage_in=FT_TRUE
   else
     OptimizeDataStorage_in=FT_FALSE
   endif
   error_out=Epetra_CrsMatrix_FillComplete_UsingMaps(this%CrsMatrix_id,DomainMap%get_EpetraMap_ID(),&
        RangeMap%get_EpetraMap_ID(),OptimizeDataStorage_in)
   if (present(err)) err=error(error_out,'Epetra_CrsMatrix%FillComplete_Map: failed.')
  end subroutine

 subroutine ExtractGlobalRowCopy(this,GlobalRow,values,indices,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: GlobalRow
   integer(c_int)                      :: NumEntries,length  ! This is just there to match the binding
   real(c_double), allocatable, dimension(:),intent(out):: values 
   integer(c_int), allocatable, dimension(:),intent(out):: indices 
   type(error), optional,       intent(out)  :: err
   integer(c_int)                            :: error_out
   allocate(values(this%NumGlobalEntries(GlobalRow)))
   allocate(indices(this%NumGlobalEntries(GlobalRow)))
   length=size(values)
   error_out=Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices(this%CrsMatrix_id,GlobalRow,&
                         length,NumEntries,values,indices)
   call assert(NumEntries==this%NumGlobalEntries(GlobalRow) &
        ,error_message('ExtractGlobalRowCopy: Mismatch between NumEntries and size(values)'))
   if (present(err)) err=error(error_out,'Epetra_CrsMatrix%ExtractGlobalRowCopy: failed')
 end subroutine

 integer(c_int) function NumMyEntries(this,MyRow)
   use iso_c_binding, only : c_int
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: MyRow
   integer(c_int)                      :: MyRow_c
   MyRow_c=MyRow-FT_Index_OffSet ! To account for Fortran index base 1
   NumMyEntries=Epetra_CrsMatrix_NumMyEntries(this%CrsMatrix_id,MyRow_c)
 end function

 integer(c_int) function NumMyRowEntries(this,MyRow)
   use iso_c_binding, only : c_int
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: MyRow
   integer(c_int)                      :: MyRow_c
   integer(c_int)                      :: junk
   MyRow_c=MyRow-FT_Index_OffSet ! To account for Fortran index base 1
   junk=Epetra_CrsMatrix_NumMyRowEntries(this%CrsMatrix_id,MyRow_c,NumMyRowEntries)
 end function

 integer(c_int) function MaxNumEntries(this)
   use iso_c_binding, only: c_int
   class(Epetra_CrsMatrix), intent(in) :: this
   MaxNumEntries=Epetra_CrsMatrix_MaxNumEntries(this%CrsMatrix_id)
 end function
 
 subroutine Multiply_Vector(this,TransA,x,y,err)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   use FEpetra_Vector, only:Epetra_Vector
   class(Epetra_CrsMatrix), intent(in) :: this
   logical, intent(in) :: TransA
   integer(FT_boolean_t)  :: TransA_in
   class(Epetra_Vector), intent(in) :: x
   class(Epetra_Vector), intent(inout) :: y 
   type(error), optional, intent(inout) :: err
   integer(c_int)                       :: error_out
   if (TransA) then
     TransA_in=FT_TRUE
   else
     TransA_in=FT_FALSE
   endif
   error_out=Epetra_CrsMatrix_Multiply_Vector(this%CrsMatrix_id,TransA_in,x%get_EpetraVector_ID(),y%get_EpetraVector_ID())    
   if (present(err)) err=error(error_out,'Epetra_CrsMatrix%Multiply_Vector: failed')
 end subroutine

 subroutine Multiply_MultiVector(this,TransA,x,y,err)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   use FEpetra_MultiVector, only:Epetra_MultiVector
   class(Epetra_CrsMatrix), intent(in) :: this
   logical, intent(in) :: TransA
   integer(FT_boolean_t)  :: TransA_in
   class(Epetra_MultiVector), intent(in) :: x
   class(Epetra_MultiVector), intent(inout) :: y 
   type(error), optional, intent(inout) :: err
   integer(c_int)                       :: error_out
   if (TransA) then
     TransA_in=FT_TRUE
   else
     TransA_in=FT_FALSE
   endif
   error_out=Epetra_CrsMatrix_Multiply_MultiVector(this%CrsMatrix_id,TransA_in,&
                                                  x%get_EpetraMultiVector_ID(),y%get_EpetraMultiVector_ID())    
   if (present(err)) err=error(error_out,'Epetra_CrsMatrix%Multiply_MultiVector: failed')
 end subroutine

 logical function MyGlobalRow(this,GID)
   use ForTrilinos_enums, only: FT_boolean_t,FT_TRUE
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int), intent(in) ::GID
   integer(FT_boolean_t) :: MyGlobalRow_out
   MyGlobalRow_out=Epetra_CrsMatrix_MyGlobalRow(this%CrsMatrix_id,GID)
   if (MyGlobalRow_out==FT_TRUE) then
     MyGLobalRow=.true.
   else
     MyGLobalRow=.false.
   endif
 end function
 
 type(Epetra_Map) function RowMatrixRowMap(this)
  use ForTrilinos_enums, only:FT_Epetra_Map_ID_t
  class(Epetra_CrsMatrix), intent(in) :: this
  type(FT_Epetra_Map_ID_t)            :: RowMap_ID
  RowMap_ID=Epetra_CrsMatrix_RowMatrixRowMap(this%CrsMatrix_id) 
  RowMatrixRowMap=Epetra_Map(RowMap_ID)
 end function

 type(Epetra_Map) function RowMap(this)
  use ForTrilinos_enums, only:FT_Epetra_Map_ID_t
  class(Epetra_CrsMatrix), intent(in) :: this
  type(FT_Epetra_Map_ID_t)            :: RowMap_ID
  RowMap_ID=Epetra_CrsMatrix_RowMap(this%CrsMatrix_id) 
  RowMap=Epetra_Map(RowMap_ID)
 end function

 integer(c_int) function NumGlobalEntries(this,GlobalRow)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int), intent(in) :: GlobalRow
   NumGlobalEntries=Epetra_CrsMatrix_NumGlobalEntries(this%CrsMatrix_id,GlobalRow)
 end function

  subroutine invalidate_EpetraCrsMatrix_ID(this)
    class(Epetra_CrsMatrix) ,intent(inout) :: this
    call this%invalidate_EpetraRowMatrix_ID
    this%CrsMatrix_id%table = FT_Invalid_ID
    this%CrsMatrix_id%index = FT_Invalid_Index
    this%CrsMatrix_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraCrsMatrix(this)
    class(Epetra_CrsMatrix) ,intent(inout) :: this
    call this%ctrilinos_delete_EpetraRowMatrix()
    call Epetra_CrsMatrix_Destroy( this%CrsMatrix_id ) 
  end subroutine
end module 

