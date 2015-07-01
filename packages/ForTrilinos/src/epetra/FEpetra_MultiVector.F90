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


module FEpetra_MultiVector
  use ForTrilinos_enums ,only: FT_Epetra_MultiVector_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_DistObject, only: Epetra_DistObject
  use FEpetra_BlockMap  ,only: Epetra_BlockMap
  use iso_c_binding     ,only: c_int,c_double,c_char
  use forepetra
  implicit none
  private                      ! Hide everything by default
  public :: Epetra_MultiVector ! Expose type/constructors/methods

  type ,extends(universal)                    :: Epetra_MultiVector 
    private
    type(FT_Epetra_MultiVector_ID_t) :: MultiVector_id 
    type(Epetra_DistObject) :: DistObject
  contains
     ! Constructors
     procedure ,private :: duplicate_
     procedure ,private :: create_
     procedure ,private :: from_struct_
     generic :: Epetra_MultiVector_ => from_struct_,duplicate_,create_
     ! Post-construction modification procedure 
     procedure ,private:: ReplaceGlobalValue_NoOffset
     procedure ,private:: ReplaceGlobalValue_BlockPos
     generic :: ReplaceGlobalValue=>ReplaceGlobalValue_NoOffset,ReplaceGlobalValue_BlockPos
     procedure ,private:: ReplaceMyValue_NoOffset
     procedure ,private:: ReplaceMyValue_BlockPos
     generic :: ReplaceMyValue=>ReplaceMyValue_NoOffset,ReplaceMyValue_BlockPos
     procedure ,private:: SumIntoGlobalValue_NoOffset
     procedure ,private:: SumIntoGlobalValue_BlockPos
     generic :: SumIntoGlobalValue=>SumIntoGlobalValue_NoOffset,SumIntoGlobalValue_BlockPos
     procedure ,private:: SumIntoMyValue_NoOffset
     procedure ,private:: SumIntoMyValue_BlockPos
     generic :: SumIntoMyValue=>SumIntoMyValue_NoOffset,SumIntoMyValue_BlockPos
     !Mathematical Methods
     procedure         :: Dot
     procedure         :: Abs
     procedure         :: Reciprocal
     procedure         :: Scale_Self
     procedure         :: Scale_Other
     generic :: Scale => Scale_Self,Scale_Other
     procedure         :: PutScalar
     procedure         :: Random
     ! Extraction Methods
     procedure         :: ExtractCopy_2DA
     procedure         :: ExtractCopy_MoldR2
     generic :: ExtractCopy=>ExtractCopy_2DA, ExtractCopy_MoldR2
     ! Mathematical Methods
     procedure         :: Update_WithA
     procedure         :: Update_WithAB
     generic :: Update => Update_WithA,Update_WithAB
     procedure         :: Norm1
     procedure         :: Norm2
     procedure         :: NormInf
     procedure         :: NormWeighted
     procedure         :: MinValue
     procedure         :: MaxValue
     procedure         :: MeanValue
     procedure         :: Multiply_Matrix
     procedure         :: Multiply_ByEl
     generic :: Multiply => Multiply_Matrix,Multiply_ByEl
     procedure         :: ReciprocalMultiply
     !Attribute Access Functions
     procedure         :: NumVectors
     procedure         :: MyLength
     procedure         :: GlobalLength
     procedure         :: Stride
     procedure         :: ConstantStride
     procedure,private :: Export_UsingImporter
     procedure,private :: Export_UsingExporter
     generic           :: export=>Export_UsingImporter,Export_UsingExporter
     procedure,private :: Import_UsingImporter
     procedure,private :: Import_UsingExporter
     generic           :: import=>Import_UsingImporter,Import_UsingExporter
     !ForTrilinos Developers only -- not to be invoked by end-user applications:
     procedure         :: invalidate_id => invalidate_EpetraMultiVector_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraMultiVector
     procedure         :: component_finalization => component_finalization_EpetraMultiVector
     procedure         :: get_EpetraMultiVector_ID 
     procedure ,nopass :: alias_EpetraMultiVector_ID
     procedure         :: generalize 
  end type

   interface Epetra_MultiVector ! constructors
     module procedure create,duplicate,from_struct
   end interface

contains

  subroutine from_struct_(this,id)
     class(Epetra_MultiVector) ,intent(out) :: this
     type(FT_Epetra_MultiVector_ID_t) ,intent(in) :: id
     this%MultiVector_id = id  
     this%DistObject = Epetra_DistObject(this%DistObject%alias_EpetraDistObject_ID(this%generalize()))
     call this%register_self
  end subroutine

  function from_struct(id) result(new_Epetra_MultiVector)
    type(Epetra_MultiVector) :: new_Epetra_MultiVector
    type(FT_Epetra_MultiVector_ID_t) ,intent(in) :: id
    call new_Epetra_MultiVector%Epetra_MultiVector_(id)
  end function

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut = true);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut ); 

  subroutine create_(this,BlockMap,Num_Vectors,zero)
    use ForTrilinos_enums ,only: FT_boolean_t,FT_TRUE,FT_FALSE
    use iso_c_binding     ,only: c_int
    class(Epetra_MultiVector) ,intent(out) :: this
    class(Epetra_BlockMap) ,intent(in) :: BlockMap
    integer(c_int)         ,intent(in) :: Num_Vectors
    logical  ,intent(in) :: zero 
    integer(FT_boolean_t) :: zero_in 
    if (zero) then
      zero_in=FT_TRUE
    else 
      zero_in=FT_FALSE
    end if
    call this%Epetra_MultiVector_(Epetra_MultiVector_Create(BlockMap%get_EpetraBlockMap_ID(),Num_Vectors,zero_in))
  end subroutine

  function create(BlockMap,Num_Vectors,zero) result(new_Epetra_MultiVector)
    use iso_c_binding     ,only: c_int
    type(Epetra_MultiVector) :: new_Epetra_MultiVector
    class(Epetra_BlockMap) ,intent(in) :: BlockMap
    integer(c_int)         ,intent(in) :: Num_Vectors
    logical  ,intent(in) :: zero 
    call new_Epetra_MultiVector%Epetra_MultiVector_(BlockMap,Num_Vectors,zero)
  end function

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_MultiVector& Source);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( CT_Epetra_MultiVector_ID_t SourceID );

  subroutine duplicate_(this,copy)
    class(Epetra_MultiVector) ,intent(in) :: this
    type(Epetra_MultiVector) ,intent(out) :: copy
    call copy%Epetra_MultiVector_(Epetra_MultiVector_Duplicate(this%MultiVector_id))
  end subroutine

  type(Epetra_MultiVector) function duplicate(original)
    type(Epetra_MultiVector) ,intent(in) :: original
    call original%Epetra_MultiVector_(duplicate)
  end function

  type(FT_Epetra_MultiVector_ID_t) function get_EpetraMultiVector_ID(this)
    class(Epetra_MultiVector) ,intent(in) :: this 
    get_EpetraMultiVector_ID=this%MultiVector_id
  end function
  
  type(FT_Epetra_MultiVector_ID_t) function alias_EpetraMultiVector_ID(generic_id)
    use ForTrilinos_table_man,only: CT_Alias
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t,FT_Epetra_MultiVector_ID
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_MultiVector_ID),stat=status)
    ierr=error(status,'FEpetra_MultiVector:alias_EpetraMultiVector_ID')
    call ierr%check_success()
    alias_EpetraMultiVector_ID=degeneralize_EpetraMultiVector(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only : c_loc
   class(Epetra_MultiVector) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%MultiVector_id))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_MultiVector) ,intent(in) ,target :: this
   ! generalize = Epetra_MultiVector_Generalize ( this%MultiVector_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_MultiVector_ID_t) function degeneralize_EpetraMultiVector(generic_id) bind(C)
  ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_MultiVector_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                      ,value   :: generic_id
    type(FT_Epetra_MultiVector_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraMultiVector = local_ptr
  ! ____ Use for ForTrilinos function implementation ______

  ! ____ Use for CTrilinos function implementation ______
  !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
  !degeneralize_EpetraMultiVector = Epetra_MultiVector_Degeneralize(generic_id)
  ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine ReplaceGlobalValue_NoOffset(this,GlobalRow,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalRow
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_c
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_c=VectorIndex-FT_Index_OffSet ! To account for Fortran index base 1
    error_out=Epetra_MultiVector_ReplaceGlobalValue(this%MultiVector_id,GlobalRow,VectorIndex_c,ScalarValue)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%ReplaceGlobalValue_NoOffset: failed')
  end subroutine
  
  subroutine ReplaceGlobalValue_BlockPos(this,GlobalBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalBlockRow
    integer(c_int), intent(in) :: BlockRowOffset
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_c ! To account for Fortran index base 1
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_c=VectorIndex-FT_Index_OffSet ! To account for Fortran index base 1
    error_out=Epetra_MultiVector_ReplaceGlobalValue_BlockPos(this%MultiVector_id,GlobalBlockRow,BlockRowOffset,&
                                                             VectorIndex_c,ScalarValue)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%ReplaceGlobalValue_BlockPos: failed')
  end subroutine
 
  subroutine ReplaceMyValue_NoOffset(this,MyRow,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyRow
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_c
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_c=VectorIndex-FT_Index_OffSet ! To account for Fortran index base 1
    error_out=Epetra_MultiVector_ReplaceMyValue(this%MultiVector_id,MyRow,VectorIndex_c,ScalarValue)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%ReplaceMyValue_NoOffset: failed')
  end subroutine
  
  subroutine ReplaceMyValue_BlockPos(this,MyBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyBlockRow
    integer(c_int), intent(in) :: BlockRowOffset
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_c ! To account for Fortran index base 1
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_c=VectorIndex-FT_Index_OffSet ! To account for Fortran index base 1
    error_out=Epetra_MultiVector_ReplaceMyValue_BlockPos(this%MultiVector_id,MyBlockRow,BlockRowOffset,&
                                                         VectorIndex_c,ScalarValue)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%ReplaceMyValue_BlockPos: failed')
  end subroutine

  subroutine SumIntoGlobalValue_NoOffset ( this,GlobalRow,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalRow
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_c ! To account for Fortran index base 1
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_c=VectorIndex-FT_Index_OffSet ! To account for Fortran index base 1
    error_out=Epetra_MultiVector_SumIntoGlobalValue(this%MultiVector_id,GlobalRow,VectorIndex_c,ScalarValue)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%SumIntoGlobalValue_NoOffset: failed')
  end subroutine
  
  subroutine SumIntoGlobalValue_BlockPos ( this,GlobalBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalBlockRow
    integer(c_int), intent(in) :: BlockRowOffset
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_c ! To account for Fortran index base 1
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_c=VectorIndex-FT_Index_OffSet ! To account for Fortran index base 1
    error_out=Epetra_MultiVector_SumIntoGlobalValue_BlockPos(this%MultiVector_id,GlobalBlockRow,BlockRowOffset,&
                                                             VectorIndex_c,ScalarValue)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%SumIntoGlobalValue_BlockPos: failed')
  end subroutine

  subroutine SumIntoMyValue_NoOffset ( this,MyRow,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyRow
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_c ! To account for Fortran index base 1
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_c=VectorIndex-FT_Index_OffSet ! To account for Fortran index base 1
    error_out=Epetra_MultiVector_SumIntoMyValue(this%MultiVector_id,MyRow,VectorIndex_c,ScalarValue)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%SumIntoMyValue_NoOffset: failed')
  end subroutine
  
  subroutine SumIntoMyValue_BlockPos ( this,MyBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyBlockRow
    integer(c_int), intent(in) :: BlockRowOffset
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_c ! To account for Fortran index base 1
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_c=VectorIndex-FT_Index_OffSet ! To account for Fortran index base 1
    error_out=Epetra_MultiVector_SumIntoMyValue_BlockPos(this%MultiVector_id,MyBlockRow,BlockRowOffset,&
                                                         VectorIndex_c,ScalarValue)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%SumIntoMyValue_BlockPos: failed')
  end subroutine


  function Dot(this,x,err) result(dot_)
    class(Epetra_MultiVector), intent(in) :: this
    class(Epetra_MultiVector), intent(in) :: x
    real(c_double),dimension(:),allocatable :: dot_ 
    type(error),optional,intent(out)   :: err
    integer(c_int)                        :: error_out
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.allocated(dot_)) then
      allocate(dot_(this%NumVectors()),stat=status)
      ierr=error(status,'FEpetra_MultiVector:Dot')
      call ierr%check_success()
    endif
    error_out=Epetra_MultiVector_Dot(this%MultiVector_id,x%MultiVector_id,dot_)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Dot: failed')
  end function 

  subroutine Abs(this,A,err)
    class(Epetra_MultiVector), intent(inout) :: this
    class(Epetra_MultiVector), intent(in) :: A 
    type(error),optional,intent(out)   :: err
    integer(c_int)                        :: error_out
    error_out=Epetra_MultiVector_Abs(this%MultiVector_id,A%MultiVector_id)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Abs: failed')
  end subroutine

  subroutine Reciprocal(this,A,err)
    class(Epetra_MultiVector), intent(inout) :: this
    class(Epetra_MultiVector), intent(in) :: A
    type(error),optional,intent(out)   :: err
    integer(c_int)                        :: error_out
    error_out=Epetra_MultiVector_Reciprocal(this%MultiVector_id,A%MultiVector_id)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Reciprocal: failed')
  end subroutine
 
  subroutine Scale_Self(this,scalar_value,err)
    class(Epetra_MultiVector), intent(inout) :: this
    real(c_double),            intent(in)    :: scalar_value
    type(error),optional,intent(out)      :: err
    integer(c_int)                           :: error_out
    error_out=Epetra_MultiVector_Scale_Self(this%MultiVector_id,scalar_value)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Scale_Self: failed')
  end subroutine
  
  subroutine Scale_Other(this,scalar_value,MultiVector,err)
    class(Epetra_MultiVector), intent(inout) :: this
    class(Epetra_MultiVector), intent(in)    :: MultiVector 
    real(c_double),            intent(in)    :: scalar_value
    type(error),optional,intent(out)      :: err
    integer(c_int)                           :: error_out
    error_out=Epetra_MultiVector_Scale(this%MultiVector_id,scalar_value,MultiVector%MultiVector_id)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Scale_Other: failed')
  end subroutine

  subroutine PutScalar(this,scalar,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    type(error) ,optional  ,intent(out)   :: err
    real(c_double)            ,intent(in)    :: scalar
    integer(c_int)                           :: error_out
    error_out=Epetra_MultiVector_PutScalar(this%MultiVector_id,scalar)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Put_Scalar: failed')
  end subroutine
 
  subroutine Random(this,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    type(error) ,optional  ,intent(out)   :: err
    integer(c_int)                           :: error_out
    error_out = Epetra_MultiVector_Random (this%MultiVector_id)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Random: failed')
  end subroutine

  function ExtractCopy_2DA(this,MyLDA,err) result(A)
    class(Epetra_MultiVector),intent(in) :: this
    real(c_double),dimension(:,:),allocatable :: A
    integer(c_int), intent(in) :: MyLDA
    type(error),optional,intent(out) :: err
    integer(c_int)           :: error_out
    integer(c_int) :: status,lda
    type(error) :: ierr

    lda = max(MyLDA,this%MyLength())
    ! Note: a function result cannot possibly be allocated at this point
    allocate(A(lda,this%NumVectors()),stat=status) ! To match a user's Fortran-style array in Trilinos
    ierr=error(status,'FEpetra_MultiVector:ExtractCopy_2DA')
    call ierr%check_success()
    error_out=Epetra_MultiVector_ExtractCopy_Fill2DA(this%MultiVector_id,A,lda)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%ExtractCopy_2DA: failed')
  end function

  function ExtractCopy_MoldR2(this,mold,err) result(A)
    class(Epetra_MultiVector), intent(in)       :: this
    real(c_double), dimension(:,:), allocatable :: A
    real(c_double), dimension(:,:), intent(in)  :: mold
    type(error),optional,intent(out) :: err
    integer(c_int) :: error_out
    integer(c_int) :: status,lda
    type(error)    :: ierr

    lda = this%MyLength()
    allocate(A(lda,this%NumVectors()),stat=status) ! To match a user's Fortran-style array in Trilinos
    ierr=error(status,'FEpetra_MultiVector:ExtractCopy_MoldR2')
    call ierr%check_success()
    error_out=Epetra_MultiVector_ExtractCopy_Fill2DA(this%MultiVector_id,A,lda)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%ExtractCopy_Fill2DA: failed')
  end function ExtractCopy_MoldR2

  subroutine Update_WithA(this,scalarA,A,scalarThis,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    real(c_double)            ,intent(in) :: scalarA
    class(Epetra_MultiVector) ,intent(in) :: A
    real(c_double)            ,intent(in) :: scalarThis
    type(error) ,optional  ,intent(out):: err
    integer(c_int)                        :: error_out
    error_out = Epetra_MultiVector_Update_WithA(this%MultiVector_id,scalarA,A%MultiVector_id,scalarThis)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Update_WithA: failed')
  end subroutine 

  subroutine Update_WithAB(this,scalarA,A,scalarB,B,scalarThis,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    real(c_double)            ,intent(in) :: scalarA
    class(Epetra_MultiVector) ,intent(in) :: A
    real(c_double)            ,intent(in) :: scalarB
    class(Epetra_MultiVector) ,intent(in) :: B
    real(c_double)            ,intent(in) :: scalarThis
    type(error) ,optional  ,intent(out):: err
    integer(c_int)                        :: error_out
    error_out = Epetra_MultiVector_Update_WithAB(this%MultiVector_id,scalarA,A%MultiVector_id,scalarB,B%MultiVector_id,scalarThis)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Update_WithAB: failed')
  end subroutine

  function Norm1(this,err) result(Norm1_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: Norm1_val 
    integer(c_int)                           :: error_out
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.allocated(Norm1_val)) then
      allocate(Norm1_val(this%NumVectors()),stat=status)
      ierr=error(status,'FEpetra_MultiVector:Norm1')
      call ierr%check_success()
    endif
    error_out = Epetra_MultiVector_Norm1(this%MultiVector_id,Norm1_val)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Norm: failed')
  end function 
 
  function Norm2(this,err) result(Norm2_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: Norm2_val 
    integer(c_int)                           :: error_out
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.allocated(Norm2_val)) then
      allocate(Norm2_val(this%NumVectors()),stat=status)
      ierr=error(status,'FEpetra_MultiVector:Norm2')
      call ierr%check_success()
    endif
    error_out = Epetra_MultiVector_Norm2(this%MultiVector_id,Norm2_val)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Norm2: failed')
  end function
 
  function NormInf(this,err) result(NormInf_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: NormInf_val 
    integer(c_int)                           :: error_out
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.allocated(NormInf_val)) then
      allocate(NormInf_val(this%NumVectors()),stat=status)
      ierr=error(status,'FEpetra_MultiVector:NormInf')
      call ierr%check_success()
    endif
    error_out = Epetra_MultiVector_NormInf(this%MultiVector_id,NormInf_val)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%NormInf: failed')
  end function 

  function NormWeighted(this,weights,err) result(NormWeighted_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    class(Epetra_MultiVector)   ,intent(in)  :: weights 
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: NormWeighted_val 
    integer(c_int)                           :: error_out
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.allocated(NormWeighted_val)) then
      allocate(NormWeighted_val(this%NumVectors()),stat=status)
      ierr=error(status,'FEpetra_MultiVector:NormWeighted')
      call ierr%check_success()
    endif
    error_out = Epetra_MultiVector_NormWeighted(this%MultiVector_id,weights%MultiVector_id,NormWeighted_val)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%NormWeighted: failed')
  end function 

  function MinValue(this,err) result(MinValue_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: MinValue_val
    integer(c_int)                           :: error_out
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.allocated(MinValue_val)) then
      allocate(MinValue_val(this%NumVectors()),stat=status)
      ierr=error(status,'FEpetra_MultiVector:MinValue')
      call ierr%check_success()
    endif
    error_out = Epetra_MultiVector_MinValue(this%MultiVector_id,MinValue_val)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%MinValue: failed')
  end function
 
 function MaxValue(this,err) result(MaxValue_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: MaxValue_val
    integer(c_int)                           :: error_out
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.allocated(MaxValue_val)) then
      allocate(MaxValue_val(this%NumVectors()),stat=status)
      ierr=error(status,'FEpetra_MultiVector:MaxValue')
      call ierr%check_success()
    endif
    error_out = Epetra_MultiVector_MaxValue(this%MultiVector_id,MaxValue_val)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%MaxValue: failed')
  end function

 function MeanValue(this,err) result(MeanValue_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: MeanValue_val
    integer(c_int)                           :: error_out
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.allocated(MeanValue_val)) then
      allocate(MeanValue_val(this%NumVectors()),stat=status)
      ierr=error(status,'FEpetra_MultiVector:MeanValue')
      call ierr%check_success()
    endif
    error_out = Epetra_MultiVector_MeanValue(this%MultiVector_id,MeanValue_val)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%MeanValue: failed')
  end function

  subroutine Multiply_Matrix(this,TransA,TransB,ScalarAB,A,B,ScalarThis,err) 
    class(Epetra_MultiVector)   ,intent(in) :: this
    character(c_char), intent(in) :: TransA
    character(c_char), intent(in) :: TransB
    class(Epetra_MultiVector)   ,intent(in)  :: A,B
    real(c_double), intent(in) :: ScalarAB, ScalarThis 
    type(error) ,optional    ,intent(out) :: err
    integer(c_int)                           :: error_out
    error_out = Epetra_MultiVector_Multiply_Matrix(this%MultiVector_id,TransA,TransB,ScalarAB,&
                                 A%MultiVector_id,B%MultiVector_id,ScalarThis)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Multiply_Matrix: failed')
  end subroutine

 subroutine Multiply_ByEl(this,ScalarAB,A,B,ScalarThis,err)
    class(Epetra_MultiVector)   ,intent(in) :: this
    class(Epetra_MultiVector)   ,intent(in)  :: A,B
    real(c_double), intent(in) :: ScalarAB, ScalarThis
    type(error) ,optional    ,intent(out) :: err
    integer(c_int)                           :: error_out
    error_out = Epetra_MultiVector_Multiply_ByEl(this%MultiVector_id,ScalarAB,&
                                 A%MultiVector_id,B%MultiVector_id,ScalarThis)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Multiply_ByEl: failed')
  end subroutine

  subroutine ReciprocalMultiply(this,ScalarAB,A,B,ScalarThis,err)
    class(Epetra_MultiVector)   ,intent(in) :: this
    class(Epetra_MultiVector)   ,intent(in)  :: A,B
    real(c_double), intent(in) :: ScalarAB, ScalarThis
    type(error) ,optional    ,intent(out) :: err
    integer(c_int)                           :: error_out
    error_out = Epetra_MultiVector_ReciprocalMultiply(this%MultiVector_id,ScalarAB,&
                                 A%MultiVector_id,B%MultiVector_id,ScalarThis)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%ReciprocalMultiply: failed')
  end subroutine

  integer(c_int) function NumVectors(this)
    class(Epetra_MultiVector) ,intent(in) :: this
    NumVectors=Epetra_MultiVector_NumVectors(this%MultiVector_id)
  end function 

  integer(c_int) function MyLength(this)
    class(Epetra_MultiVector) ,intent(in) :: this
    MyLength=Epetra_MultiVector_MyLength(this%MultiVector_id)
  end function 

  integer(c_int) function GlobalLength(this)
    class(Epetra_MultiVector) ,intent(in) :: this
    GlobalLength=Epetra_MultiVector_GlobalLength(this%MultiVector_id)
  end function 

  integer(c_int) function Stride(this)
    class(Epetra_MultiVector) ,intent(in) :: this
    Stride=Epetra_MultiVector_Stride(this%MultiVector_id)
  end function 

  function ConstantStride(this) result(is_constant)
    use ForTrilinos_enums ,only: FT_boolean_t
    integer(FT_boolean_t) :: is_constant
    class(Epetra_MultiVector) ,intent(in) :: this
    is_constant=Epetra_MultiVector_ConstantStride(this%MultiVector_id)
  end function 

 subroutine Export_UsingExporter(this,A,exporter,CombineMode,indexor,err)
    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
    use FEpetra_Export, only: Epetra_Export
    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
    class(Epetra_MultiVector), intent(in) :: this,A 
    type(Epetra_Export),intent(in) :: exporter
    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
    type(Epetra_OffsetIndex), intent(in) :: indexor
    type(error),optional,intent(out) :: err 
    integer(c_int)     :: error_out
    call this%DistObject%export(A%DistObject,exporter,CombineMode,indexor,err)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Export_UsingExporter: failed')
  end subroutine

  subroutine Export_UsingImporter(this,A,importer,CombineMode,indexor,err)
    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
    use FEpetra_Import, only: Epetra_Import
    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
    class(Epetra_MultiVector), intent(in) :: this,A
    type(Epetra_Import),intent(in) :: importer
    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
    type(Epetra_OffsetIndex), intent(in) :: indexor
    type(error),optional,intent(out) :: err 
    integer(c_int)     :: error_out
    call this%DistObject%export(A%DistObject,importer,CombineMode,indexor,err)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Export_UsingImporter: failed')
  end subroutine

 subroutine Import_UsingExporter(this,A,exporter,CombineMode,indexor,err)
    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
    use FEpetra_Export, only: Epetra_Export
    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
    class(Epetra_MultiVector), intent(in) :: this,A 
    type(Epetra_Export),intent(in) :: exporter
    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
    type(Epetra_OffsetIndex), intent(in) :: indexor
    type(error),optional,intent(out) :: err 
    integer(c_int)     :: error_out
    call this%DistObject%import(A%DistObject,exporter,CombineMode,indexor,err)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Import_UsingExporter: failed')
  end subroutine

  subroutine Import_UsingImporter(this,A,importer,CombineMode,indexor,err)
    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
    use FEpetra_Import, only: Epetra_Import
    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
    class(Epetra_MultiVector), intent(in) :: this,A
    type(Epetra_Import),intent(in) :: importer
    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
    type(Epetra_OffsetIndex), intent(in) :: indexor
    type(error),optional,intent(out) :: err 
    integer(c_int)     :: error_out
    call this%DistObject%import(A%DistObject,importer,CombineMode,indexor,err)
    if (present(err)) err=error(error_out,'Epetra_MultiVector%Import_UsingImporter: failed')
  end subroutine

  subroutine invalidate_EpetraMultiVector_ID(this)
    class(Epetra_MultiVector),intent(inout) :: this
    !call this%DistObject%invalidate_id
    this%MultiVector_id%table = FT_Invalid_ID
    this%MultiVector_id%index = FT_Invalid_Index 
    this%MultiVector_id%is_const = FT_FALSE
  end subroutine

  subroutine component_finalization_EpetraMultiVector(this)
    class(Epetra_MultiVector),intent(inout) :: this
    call this%DistObject%force_finalize() 
  end subroutine

  subroutine ctrilinos_delete_EpetraMultiVector(this)
    class(Epetra_MultiVector),intent(inout) :: this
    call Epetra_MultiVector_Destroy( this%MultiVector_id ) 
  end subroutine

end module 

