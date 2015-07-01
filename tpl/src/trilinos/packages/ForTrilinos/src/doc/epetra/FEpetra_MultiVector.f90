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
  use ForTrilinos_enums ! ,only: FT_Epetra_MultiVector_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
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

  type  Epetra_MultiVector !,extends(universal) "shell"
  contains
  end type

contains

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut = true);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut ); 

  !> @name Constructor Functions
  !! @{
 
  !> <BR> Epetra_MultiVector constructor conformal to a BlockMap, optionally zero the newly created vector.
  type(Epetra_MultiVector) function Epetra_MultiVector(BlockMap,Num_Vectors,zero)
    class(Epetra_BlockMap) ,intent(in) :: BlockMap &
         !< In The map to which the vector will conform 
    integer(c_int)         ,intent(in) :: Num_Vectors &
    !< In  Number of vectors in multivector
    logical  ,intent(in) :: zero  &
    !< In Optionally zero out the output. 
  end function

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_MultiVector& Source);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( CT_Epetra_MultiVector_ID_t SourceID );

  !> @name Constructor Functions
  !! @{
 
  !> <BR> Epetra_MultiVector copy constructor.
  type(Epetra_MultiVector) function Epetra_MultiVector(this)
    class(Epetra_MultiVector) ,intent(in) :: this
  end function 
  

  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces value at location (GlobalRow,VectorIndex) with ScalarValue
  !> @brief The index of the specified location must correspond to an index owned by the map on the calling process; i.e. no communication takes place.
  subroutine ReplaceGlobalValue(this,GlobalRow,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalRow &
         !< In The global row to be set 
    integer(c_int), intent(in) :: VectorIndex &
    !< In The vector index within the multivector to be set
    real(c_double), intent(in) :: ScalarValue &
    !< In The value to be set
    type(error),optional,intent(out) :: err &
         !< Returns  error information.

  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Replaces value at location (GlobalBlockRow,BlockRowOffset,VectorIndex) with ScalarValue
  !> @brief The index of the specified location must correspond to an index owned by the map on the calling process; i.e. no communication takes place.
  subroutine ReplaceGlobalValue(this,GlobalBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalBlockRow &
         !< In The global block row to be set
    integer(c_int), intent(in) :: BlockRowOffset &
         !< In The global block row offest to be set
    integer(c_int), intent(in) :: VectorIndex  &
         !< In The vector index within the multivector to be set
    real(c_double), intent(in) :: ScalarValue  &
         !< In The scalar value to be set 
    type(error),optional,intent(out) :: err &
         !< Returns  error information.
  end subroutine
 
  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Replaces value at location (MyRow,VectorIndex) with ScalarValue
  !> @brief The index of the specified location must be that of a locally owned element. 
  subroutine ReplaceMyValue(this,MyRow,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyRow &
         !< In The local row to be set
    integer(c_int), intent(in) :: VectorIndex   &
         !< In The vector index within the multivector to be set
    real(c_double), intent(in) :: ScalarValue  &
         !< In The scalar value to be set  
    type(error),optional,intent(out) :: err &
         !< Returns  error information.
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces value at location (MyBlockRow,BlockRowOffset,VectorIndex) with ScalarValue
  subroutine ReplaceMyValue(this,MyBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyBlockRow &
         !< In The local block row to be set 
    integer(c_int), intent(in) :: BlockRowOffset &
         !< In The local block row offset to be set 
    integer(c_int), intent(in) :: VectorIndex   &
         !< In The vector index within the multivector to be set
    real(c_double), intent(in) :: ScalarValue  &
         !< In The scalar value to be set  
    type(error),optional,intent(out) :: err
  end subroutine


  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Adds ScalarValue to value at location (GlobalRow,VectorIndex) 
  subroutine SumIntoGlobalValue ( this,GlobalRow,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalRow &
         !< In The global row to be modified
    integer(c_int), intent(in) :: VectorIndex   &
         !< In The vector index within the multivector to be modified
    real(c_double), intent(in) :: ScalarValue   &
         !< In The scalar value to be added 
    type(error),optional,intent(out) :: err
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Adds ScalarValue to value at location (GlobalBlockRow,BlockRowOffset,VectorIndex) 
  subroutine SumIntoGlobalValue ( this,GlobalBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalBlockRow &
         !< In The global block row to be modified
    integer(c_int), intent(in) :: BlockRowOffset &
         !< In The global block row offset to be modified
    integer(c_int), intent(in) :: VectorIndex  &
         !< In The vector index within the multivector to be modified
    real(c_double), intent(in) :: ScalarValue   &
         !< In The scalar value to be added 
    type(error),optional,intent(out) :: err &
         !< Returns  error information.
  end subroutine

  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Adds ScalarValue to value at location (MyRow,VectorIndex) 
  subroutine SumIntoMyValue ( this,MyRow,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyRow &
         !< In The local row to be modified
    integer(c_int), intent(in) :: VectorIndex  &
         !< In The vector index within the multivector to be modified
    real(c_double), intent(in) :: ScalarValue   &
         !< In The scalar value to be added  
    type(error),optional,intent(out) :: err &
         !< Returns  error information.
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Adds ScalarValue to value at location (MyBlockRow,BlockRowOffset,VectorIndex) 
  subroutine SumIntoMyValue ( this,MyBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyBlockRow &
         !< In The local block row to be modified
    integer(c_int), intent(in) :: BlockRowOffset &
         !< In The local block row offset to be modified
    integer(c_int), intent(in) :: VectorIndex  &
         !< In The vector index within the multivector to be modified
    real(c_double), intent(in) :: ScalarValue  &
         !< In The scalar value to be added  
    type(error),optional,intent(out) :: err &
         !< Returns  error information.
  end subroutine


  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Replaces all entries with scalar value.
  subroutine PutScalar(this,scalar,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    real(c_double)            ,intent(in)    :: scalar &
         !< In The scalar to which all entries will be set
    type(error) ,optional  ,intent(out)   :: err &
         !< Returns  error information.
  end subroutine
 
  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces all entries with random values. 
  subroutine Random(this,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    type(error) ,optional  ,intent(out)   :: err &
         !< Returns  error information.
  end subroutine



  !> @name Mathematical methods
  !! @{

  !> <BR> Computes the scalar product of corresponding pairs of vectors.
  function Dot(this,x,err) result(dot_)
    class(Epetra_MultiVector), intent(in) :: this
    class(Epetra_MultiVector), intent(in) :: x &
         !< In The multivector to be used in conjunction with the "this"  multivector
    real(c_double),dimension(:),allocatable :: dot_  &
         !< Result: dot_(i) will contain the dot product of corresponding i-th vectors. 
    type(error),optional,intent(out)   :: err &
         !< Returns  error information.
  end function 

  !> @name Mathematical methods
  !! @{

  !> <BR> Replaces target with element-wise absolute value of input
  subroutine Abs(this,A,err)
    class(Epetra_MultiVector), intent(inout) :: this
    class(Epetra_MultiVector), intent(in) :: A &
         !< In The source multivector.
    type(error),optional,intent(out)   :: err &
         !< Returns  error information.
  end subroutine

  !> @name Mathematical methods
  !! @{

  !> <BR> Reciprocal  replaces target with element-wise reciprocal value of input
  subroutine Reciprocal(this,A,err)
    class(Epetra_MultiVector), intent(inout) :: this
    class(Epetra_MultiVector), intent(in) :: A &
         !< In The source multivector.
    type(error),optional,intent(out)   :: err &
         !< Returns  error information.
  end subroutine
 
  !> @name Mathematical methods
  !! @{

  !> <BR> Scales current values  this = scalar_value*this
  subroutine Scale(this,scalar_value,err)
    class(Epetra_MultiVector), intent(inout) :: this
    real(c_double),            intent(in)    :: scalar_value &
         !< In The scale factor.
    type(error),optional,intent(out)      :: err &
         !< Returns  error information.
  end subroutine
  
  !> @name Mathematical methods
  !! @{

  !> <BR> Replaces current values with scaled input  this = scalar_value*MultiVector
  subroutine Scale(this,scalar_value,MultiVector,err)
    class(Epetra_MultiVector), intent(inout) :: this
    class(Epetra_MultiVector), intent(in)    :: MultiVector  &
         !< In The source multivector.
    real(c_double),            intent(in)    :: scalar_value &
         !< In The scale factor.
    type(error),optional,intent(out)      :: err &
         !< Returns  error information.
  end subroutine


  !> @name Mathematical methods
  !! @{

  !> <BR> Computes 1-norm of each vector in the input
  function Norm1(this,err) result(Norm1_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    real(c_double) ,dimension(:),allocatable :: Norm1_val  &
         !< Result: Norm1_val(i) will contain the 1-norm of the i-th vector. 
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end function 
 
  !> @name Mathematical methods
  !! @{

  !> <BR> Computes 2-norm of each vector in the input
  function Norm2(this,err) result(Norm2_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    real(c_double) ,dimension(:),allocatable :: Norm2_val  &
         !< Result: Norm2_val(i) will contain the 2-norm of the i-th vector. 
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end function
 
  !> @name Mathematical methods
  !! @{

  !> <BR> Computes infinity norm of each vector in the input
  function NormInf(this,err) result(NormInf_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    real(c_double) ,dimension(:),allocatable :: NormInf_val  &
         !< Result: NormInf_val(i) will contain the infinity-norm of the i-th vector. 
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end function 

  !> @name Mathematical methods
  !! @{

  !> <BR> Computes weighted  norm (RMS norm) of each vector in the input
  function NormWeighted(this,weights,err) result(NormWeighted_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    class(Epetra_MultiVector)   ,intent(in)  :: weights   &
         !< In: the weights.
    real(c_double) ,dimension(:),allocatable :: NormWeighted_val   &
         !< Result: NormWieghted_val(i) will contain the RMS-norm of the i-th vector. 
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end function 

  !> @name Mathematical methods
  !! @{

  !> <BR> Computes minimum value of each vector in the input
  function MinValue(this,err) result(MinValue_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    real(c_double) ,dimension(:),allocatable :: MinValue_val  &
         !< Result: MinValue_val(i) will contain the minimum value of the i-th vector. 
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end function
 
  !> @name Mathematical methods
  !! @{

  !> <BR> MaxValue: compute maximum value of each vector in the input
 function MaxValue(this,err) result(MaxValue_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    real(c_double) ,dimension(:),allocatable :: MaxValue_val &
         !< Result: MaxValue_val(i) will contain the maximum value of the i-th vector. 
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end function

  !> @name Mathematical methods
  !! @{

  !> <BR> Computes mean (average) value of each vector in the input
 function MeanValue(this,err) result(MeanValue_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    real(c_double) ,dimension(:),allocatable :: MeanValue_val &
         !< Result: MeanValue_val(i) will contain the average value of the i-th vector. 
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end function

  !> @name Mathematical methods
  !! @{

  !> <BR> Matrix-matrix multiplication  this = ScalarThis*This + ScalarAB*MATMUL(A,B)
    !> @brief This routine will compute the product of the two multivectors A and B and use it to update "this", like the BLAS routine GEMM. 
  subroutine Multiply(this,TransA,TransB,ScalarAB,A,B,ScalarThis,err) 
    class(Epetra_MultiVector)   ,intent(in) :: this
    character(c_char), intent(in) :: TransA &
    !< In: Choose transpose status of A
    character(c_char), intent(in) :: TransB &
    !< In: Choose transpose status of B
    class(Epetra_MultiVector)   ,intent(in)  :: A &
    !< In: Input multivector A 
    class(Epetra_MultiVector)   ,intent(in)  :: B &
    !< In: Input multivector B
    real(c_double), intent(in) :: ScalarAB &
    !< In: scale factor for product MATMUL(A,B)
    real(c_double), intent(in) :: ScalarThis  &
    !< In: scale factor for "this"
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end subroutine

  !> @name Mathematical methods
  !! @{

  !> <BR> Element-by-element multiplication  this = ScalarThis*This + ScalarAB*A*B
  subroutine Multiply(this,ScalarAB,A,B,ScalarThis,err)
    class(Epetra_MultiVector)   ,intent(in) :: this
    class(Epetra_MultiVector)   ,intent(in)  :: A &
    !< In: Input multivector A 
    class(Epetra_MultiVector)   ,intent(in)  :: B &
    !< In: Input multivector B
    real(c_double), intent(in) :: ScalarAB &
    !< In: scale factor for product A*B
    real(c_double), intent(in) :: ScalarThis  &
    !< In: scale factor for "this"
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end subroutine

  !> @name Mathematical methods
  !! @{

  !> <BR> Element-by-element multiplication by reciprocal   this = ScalarThis*This + ScalarAB*A/B
  subroutine ReciprocalMultiply(this,ScalarAB,A,B,ScalarThis,err)
    class(Epetra_MultiVector)   ,intent(in) :: this
    class(Epetra_MultiVector)   ,intent(in)  :: A &
    !< In: Input multivector A 
    class(Epetra_MultiVector)   ,intent(in)  :: B &
    !< In: Input multivector B
    real(c_double), intent(in) :: ScalarAB &
    !< In: scale factor for reciprocal product A/B
    real(c_double), intent(in) :: ScalarThis  &
    !< In: scale factor for "this"
    type(error) ,optional    ,intent(out) :: err &
         !< Returns  error information.
  end subroutine


  !> @name Mathematical methods
  !! @{

  !> <BR> Updates with scaled copy of input this = ScalarThis*This + ScalarA*A
  subroutine Update(this,scalarA,A,scalarThis,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    class(Epetra_MultiVector)   ,intent(in)  :: A &
    !< In: Input multivector A 
    real(c_double), intent(in) :: ScalarA &
    !< In: scale factor for p A
    real(c_double), intent(in) :: ScalarThis  &
    !< In: scale factor for "this"
    type(error) ,optional  ,intent(out):: err &
         !< Returns  error information.
  end subroutine 

  !> @name Mathematical methods
  !! @{

  !> <BR> Updates with scaled copies of inputs this = ScalarThis*This + ScalarA*A + ScalarB*B
  subroutine Update(this,scalarA,A,scalarB,B,scalarThis,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    real(c_double), intent(in) :: ScalarA &
    !< In: scale factor for A
    class(Epetra_MultiVector)   ,intent(in)  :: A &
    !< In: Input multivector A 
    real(c_double), intent(in) :: ScalarB &
    !< In: scale factor for B
    class(Epetra_MultiVector)   ,intent(in)  :: B &
    !< In: Input multivector B
    real(c_double), intent(in) :: ScalarThis  &
    !< In: scale factor for "this"
    type(error) ,optional  ,intent(out):: err &
         !< Returns  error information.
  end subroutine


  !> @name Extraction methods
  !! @{

  !> <BR> Copies multivector contents into target
  !> @brief The input argument MyLDA is a user request for the size of the output; the actual size will be the maximum between this and the stride of the multivector object.
  function ExtractCopy(this,MyLDA,err) result(A)
    class(Epetra_MultiVector),intent(in) :: this
    integer(c_int),intent(in) :: MyLDA &
         !< In: Minimum leading dimension of result. 
    real(c_double),dimension(:,:),allocatable :: A &
         !< Result: the copy of the local data
    type(error),optional,intent(out) :: err &
         !< Returns  error information.
  end function


  !> @name Attribute access
  !! @{

  !> <BR> Number of vectors in  multivector
  integer(c_int) function NumVectors(this)
    class(Epetra_MultiVector) ,intent(in) :: this
  end function 

  !> @name Attribute access
  !! @{

  !> <BR> Local vector length
  integer(c_int) function MyLength(this)
    class(Epetra_MultiVector) ,intent(in) :: this
  end function 

  !> @name Attribute access
  !! @{

  !> <BR> Global vector length
  integer(c_int) function GlobalLength(this)
    class(Epetra_MultiVector) ,intent(in) :: this
  end function 

  !> @name Attribute access
  !! @{

  !> <BR> Stride between successive vectors in multivector (only meaningful if ConstantStride()==true)
  integer(c_int) function Stride(this)
    class(Epetra_MultiVector) ,intent(in) :: this
  end function 

  !> @name Attribute access
  !! @{

  !> <BR> True if stride between successive vectors is constant
  integer(FT_boolean_t) function ConstantStride(this)
    use ForTrilinos_enums ,only:FT_Epetra_MultiVector_ID_t,FT_boolean_t
    class(Epetra_MultiVector) ,intent(in) :: this
  end function 

! !$ subroutine Export_UsingExporter(this,A,exporter,CombineMode,indexor,err)
! !$    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
! !$    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
! !$    use FEpetra_Export, only: Epetra_Export
! !$    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
! !$    type(Epetra_MultiVector), intent(in) :: this,A 
! !$    type(Epetra_Export),intent(in) :: exporter
! !$    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
! !$    type(Epetra_OffsetIndex), intent(in) :: indexor
! !$    type(error),optional,intent(out) :: err 
! !$    integer(c_int)     :: error_out
! !$    call this%DistObject%export(A%DistObject,exporter,CombineMode,indexor,err)
! !$    if (present(err)) err=error(error_out,'Epetra_MultiVector%Export_UsingExporter: failed')
! !$  end subroutine
! !$
! !$  subroutine Export_UsingImporter(this,A,importer,CombineMode,indexor,err)
! !$    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
! !$    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
! !$    use FEpetra_Import, only: Epetra_Import
! !$    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
! !$    type(Epetra_MultiVector), intent(in) :: this,A
! !$    type(Epetra_Import),intent(in) :: importer
! !$    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
! !$    type(Epetra_OffsetIndex), intent(in) :: indexor
! !$    type(error),optional,intent(out) :: err 
! !$    integer(c_int)     :: error_out
! !$    call this%DistObject%export(A%DistObject,importer,CombineMode,indexor,err)
! !$    if (present(err)) err=error(error_out,'Epetra_MultiVector%Export_UsingImporter: failed')
! !$  end subroutine
! !$
! !$ subroutine Import_UsingExporter(this,A,exporter,CombineMode,indexor,err)
! !$    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
! !$    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
! !$    use FEpetra_Export, only: Epetra_Export
! !$    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
! !$    type(Epetra_MultiVector), intent(in) :: this,A 
! !$    type(Epetra_Export),intent(in) :: exporter
! !$    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
! !$    type(Epetra_OffsetIndex), intent(in) :: indexor
! !$    type(error),optional,intent(out) :: err 
! !$    integer(c_int)     :: error_out
! !$    call this%DistObject%import(A%DistObject,exporter,CombineMode,indexor,err)
! !$    if (present(err)) err=error(error_out,'Epetra_MultiVector%Import_UsingExporter: failed')
! !$  end subroutine
! !$
! !$  subroutine Import_UsingImporter(this,A,importer,CombineMode,indexor,err)
! !$    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
! !$    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
! !$    use FEpetra_Import, only: Epetra_Import
! !$    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
! !$    type(Epetra_MultiVector), intent(in) :: this,A
! !$    type(Epetra_Import),intent(in) :: importer
! !$    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
! !$    type(Epetra_OffsetIndex), intent(in) :: indexor
! !$    type(error),optional,intent(out) :: err 
! !$    integer(c_int)     :: error_out
! !$    call this%DistObject%import(A%DistObject,importer,CombineMode,indexor,err)
! !$    if (present(err)) err=error(error_out,'Epetra_MultiVector%Import_UsingImporter: failed')
! !$  end subroutine
! !$

end module 

