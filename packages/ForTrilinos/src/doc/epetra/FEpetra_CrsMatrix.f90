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
  use iso_c_binding     ,only : c_int,c_long,c_double,c_char
  use forepetra
  implicit none
  !private                         ! Hide everything by default
  !public :: Epetra_CrsMatrix ! Expose type/constructors/methods

  !> <BR> Epetra_CrsMatrix: A class for constructing and using real-valued double-precision sparse compressed row matrices.
  !> @{

  !> @brief The Epetra_CrsMatrix class is a sparse compressed row matrix object. This matrix can be used in a parallel setting, with data distribution described by Epetra_Map attributes. The structure or graph of the matrix is defined by an Epetra_CrsGraph attribute.
   !! In addition to coefficient access, the primary operations provided by Epetra_CrsMatrix are matrix times vector and matrix times multi-vector multiplication.
   !! \n
   !! Epetra_CrsMatrix matrices can be square or rectangular.


  type ,extends(Epetra_RowMatrix)  :: Epetra_CrsMatrix !"shell"
  contains
    !Insertion/Replace/SumInto methods
     !procedure         :: PutScalar
     !procedure         :: Scale
     !procedure         :: InsertGlobalValues
     !procedure         :: ReplaceGlobalValues
    !Transformation methods
     !procedure,private :: FillComplete_Op
     !procedure,private :: FillComplete_Map
     !generic           :: FillComplete=>FillComplete_Op, FillComplete_Map 
     !Matrix data extraction routines
     !procedure         :: ExtractGlobalRowCopy
     !procedure         :: NumMyRowEntries
     !procedure         :: NumMyEntries
     !procedure         :: MaxNumEntries
    !Computational Methods
     !procedure         :: Multiply_Vector
     !procedure         :: Multiply => Multiply_MultiVector
     !Attribute Accessor Methods
     !procedure         :: RowMatrixRowMap 
     !procedure         :: RowMap
     !procedure         :: NumGlobalEntries
     !Local/Global ID method
     !procedure         :: MyGlobalRow
  end type

contains

  ! Original C++ prototype:
  ! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const int* NumEntriesPerRow,
  ! bool StaticProfile = false);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( CT_Epetra_DataAccess_E_t CV,
  !CT_Epetra_Map_ID_t RowMapID, const int * NumEntriesPerRow, boolean StaticProfile);

  !> @name Constructor Functions
  !> @{

  !> @brief Epetra_CrsMatrix constructor with variable number of indices per row.
  !! Creates a Epetra_CrsMatrix object and allocates storage.
  type(Epetra_CrsMatrix) function Epetra_CrsMatrix(CV,Row_Map,NumEntriesPerRow,StaticProfile)
    use ForTrilinos_enums,         only: FT_boolean_t,FT_FALSE,FT_TRUE
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV &
    !< (In) An FT_Epetra_DataAccess_E_t enumerated type set to Copy or View.
    class(Epetra_Map),              intent(in) :: Row_Map &
    !< (In) An Epetra_Map defining the numbering and distribution of matrix rows.
    integer(c_int), dimension(:),   intent(in) :: NumEntriesPerRow  &
    !< (In) An integer array of length NumRows such that NumEntriesPerRow(i) indicates the (approximate if StaticProfile=false) number of entries in the ith row.
    logical,        optional                   :: StaticProfile &
    !<  (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernels.
  end function

  ! Original C++ prototype:
  ! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int NumEntriesPerRow, bool StaticProfile = false);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  ! int NumEntriesPerRow, boolean StaticProfile );

  !> @name Constructor Functions
  !> @{

  !> @brief Epetra_CrsMatrix constructor with fixed number of indices per row.
  !! Creates a Epetra_CrsMatrix object and allocates storage.
  type(Epetra_CrsMatrix) function Epetra_CrsMatrix(CV,Row_Map,NumEntriesPerRow,StaticProfile)
    use ForTrilinos_enums,         only: FT_boolean_t,FT_FALSE,FT_TRUE
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV &
    !< (In) An FT_Epetra_DataAccess_E_t enumerated type set to Copy or View.
    class(Epetra_Map),              intent(in) :: Row_Map &
    !< (In) An Epetra_Map defining the numbering and distribution of matrix rows.
    integer(c_int),                 intent(in) :: NumEntriesPerRow &
    !< (In) An integer that indicates the (approximate) number of entries in the each row. Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
    logical,        optional                   :: StaticProfile &
    !< (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernels.
  end function

  ! Original C++ prototype:
  ! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, const int* NumEntriesPerRow, 
  ! bool StaticProfile = false);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( CT_Epetra_DataAccess_E_t CV, 
  !            CT_Epetra_Map_ID_t RowMapID, CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, boolean StaticProfile );

  !> @name Constructor Functions
  !> @{

  !> @brief Epetra_CrsMatrix constructor with variable number of indices per row.
  !! Creates a Epetra_CrsMatrix object and allocates storage.
  type(Epetra_CrsMatrix) function Epetra_CrsMatrix(CV,Row_Map,Col_Map,NumEntriesPerRow,StaticProfile)
    use ForTrilinos_enums,         only: FT_boolean_t,FT_FALSE,FT_TRUE
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV &
    !< (In) An FT_Epetra_DataAccess_E_t enumerated type set to Copy or View.
    class(Epetra_Map),              intent(in) :: Row_Map &
    !< (In) An Epetra_Map defining the numbering and distribution of matrix rows.
    class(Epetra_Map),              intent(in) :: Col_Map &
    !< (In) An Epetra_Map defining the set of column-indices that appear in each processor's locally owned matrix rows.
    integer(c_int), dimension(:) ,  intent(in) :: NumEntriesPerRow &
    !< (In) An integer array of length NumRows such that NumEntriesPerRow(i) indicates the (approximate if StaticProfile=false) number of entries in the ith row.
    logical,        optional                   :: StaticProfile &
    !< (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernels.
  end function

  !Original C++ prototype:
  ! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, int NumEntriesPerRow, 
  !                  bool StaticProfile = false);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  !                  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, boolean StaticProfile );

  !> @name Constructor Functions
  !> @{

  !> @brief Epetra_CrsMatrix constructor with fixed number of indices per row.
  !! Creates a Epetra_CrsMatrix object and allocates storage.
  type(Epetra_CrsMatrix) function Epetra_CrsMatrix(CV,Row_Map,Col_Map,NumEntriesPerRow,StaticProfile)
    use ForTrilinos_enums,         only: FT_boolean_t,FT_FALSE,FT_TRUE
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV &
    !< (In) An FT_Epetra_DataAccess_E_t enumerated type set to Copy or View.
    class(Epetra_Map),              intent(in) :: Row_Map &
    !< (In) An Epetra_Map defining the numbering and distribution of matrix rows.
    class(Epetra_Map),              intent(in) :: Col_Map &
    !< (In) An Epetra_Map defining the set of column-indices that appear in each processor's locally owned matrix rows.
    integer(c_int),                 intent(in) :: NumEntriesPerRow &
    !< (In) An integer that indicates the (approximate if StaticProfile=false) number of entries in the each row. Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
    logical,        optional                   :: StaticProfile &
    !< (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernels.
  end function

  
  ! Original C++ prototype:
  ! Epetra_CrsMatrix(const Epetra_CrsMatrix& RowMatrix);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( CT_Epetra_BasicRowMatrix_ID_t RowMatrixID );

  !> @name Constructor Functions
  !> @{

  !> @brief Epetra_CrsMatrix copy constructor.
  !! Creates a copy of a Epetra_CrsMatrix object.
  type(Epetra_CrsMatrix) function Epetra_CrsMatrix(this)
    type(Epetra_CrsMatrix) ,intent(in) :: this 
  end function

  !> @name Insertion/Replace/SumInto methods
  !> @{

  !> @brief Initialize all values in the matrix with constant value.
  subroutine PutScalar(this,scalar,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   real(c_double)           ,intent(in):: scalar &
   !< (In) Value to use.
   type(error), optional, intent(out) :: err &
   !< (Out) Integer error code, set to 0 if successful.
  end subroutine
  
  !> @name Insertion/Replace/SumInto methods
  !> @{

  !> @brief Multiply all values in the matrix by a constant value (in place: A <- ScalarConstant * A).
  subroutine Scale(this,ScalarConstant,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   real(c_double)           ,intent(in):: ScalarConstant &
   !< (In) Value to use.
   type(error), optional, intent(out) :: err &
   !< (Out) Integer error code, set to 0 if successful.
  end subroutine

  !> @name Insertion/Replace/SumInto methods
  !> @{

  !> @brief Insert a list of elements in a given global row of the matrix.
  !! This method is used to construct a matrix for the first time.  It cannot be used if the matrix structure has already been fixed (via a call to FillComplete()).
  subroutine InsertGlobalValues(this,GlobalRow,values,indices,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: GlobalRow &
   !< (In) Row number (in global coordinates) to put elements.
   real(c_double),dimension(:),intent(in):: values  &
   !< Values to enter.
   integer(c_int),dimension(:),intent(in):: indices  &
   !< (In) Global column indices corresponding to values.
   type(error), optional, intent(out) :: err &
   !< (Out) Integer error code, set to 0 if successful.
  end subroutine

  !> @name Insertion/Replace/SumInto methods
  !> @{

  !> @brief Replace specified existing values with this list of entries for a given global row of the matrix.
  subroutine ReplaceGlobalValues(this,GlobalRow,values,indices,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: GlobalRow &
   !< (In) Row number (in global coordinates) to put elements.
   real(c_double), dimension(:)        :: values  &
   !< (In) Values to enter.
   integer(c_int),    dimension(:)        :: indices  &
   !< (In) Global column indices corresponding to values.
   type(error), optional, intent(out) :: err &
   !<  Integer error code, set to 0 if successful.
  end subroutine
 
  !> @name Transformation methods
  !> @{

  !> @brief Signal that data entry is complete.  Perform transformations to local index space.
  !! This version of FillComplete assumes that the domain and range distributions are identical to the matrix row distributions. 
  subroutine FillComplete(this,OptimizeDataStorage,err)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   class(Epetra_CrsMatrix), intent(in) :: this
   logical,optional,        intent(in) :: OptimizeDataStorage &
   !< (In) If true, storage will be packed for optimal performance. 
   type(error),optional,intent(out)    :: err &
   !< error code, 0 if successful. Returns a positive warning code of 3 if the matrix is rectangular (meaning that the other overloading of FillComplete should have been called, with differen domain-map and range-map specified).
  end subroutine

  !> @name Transformation methods
  !> @{

  !> @brief Signal that data entry is complete.  Perform transformations to local index space.
  !! This version of FillComplete requires the explicit specification of the domain and range distribution maps.  These maps are used for importing and exporting vector and multi-vector elements that are needed for distributed matrix computations.
  subroutine FillComplete(this,DomainMap,RangeMap,OptimizeDataStorage,err)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   class(Epetra_CrsMatrix), intent(in) :: this
   class(Epetra_Map) , intent(in) :: DomainMap &
   !< (In) Map that describes the distribution of vector and multi-vectors in the matrix domain.
   class(Epetra_Map) , intent(in) :: RangeMap &
   !< (In) Map that describes the distribution of vector and multi-vectors in the matrix range.
   logical,optional,        intent(in) :: OptimizeDataStorage &
   !< (In) If true, storage will be packed for optimal performance.By default storage will be optimized.  If you cannot tolerate the increased temporary memory use, should set this value to false.
   type(error),optional,intent(out)    :: err &
   !< (Out) error code, 0 if successful. positive warning code of 2 if it is detected that the matrix-graph got out of sync since this matrix was constructed (for instance if graph. FillComplete() was called by another matrix that shares the graph)
  end subroutine

  !> @name Extraction Methods
  !> @{

  !> @brief Returns a copy of the specified global row in user-provided arrays.
 subroutine ExtractGlobalRowCopy(this,GlobalRow,values,indices,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: GlobalRow & 
   !< (In) Global row to extract.
   real(c_double), allocatable, dimension(:),intent(out):: values &
   !< (Out) Extracted values for this row.
   integer(c_int), allocatable, dimension(:),intent(out):: indices &
   !<  (Out) Extracted global column indices for the corresponding values.
   type(error), optional,       intent(out)  :: err &
   !< Integer error code, set to 0 if successful, non-zero if global row is not owned by calling process or if the number of entries in this row exceed the Length parameter.
 end subroutine

 !> @name Extraction Methods
 !> @{

 !> @brief Returns the current number of nonzero entries in specified local row on this processor.
 integer(c_int) function NumMyEntries(this,MyRow)
   use iso_c_binding, only : c_int
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: MyRow
 end function

 !> @name Additional methods required to implement Epetra_RowMatrix interface
 !> @{

 !> @brief Return the current number of values stored for the specified local row.
 integer(c_int) function NumMyRowEntries(this,MyRow)
   use iso_c_binding, only : c_int
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: MyRow &
   !< (In) Local row.
 end function

 !> @name Extraction Methods
 !> @{

 !> @brief Returns the maximum number of nonzero entries across all rows on this processor.
 integer(c_int) function MaxNumEntries(this)
   use iso_c_binding, only: c_int
   class(Epetra_CrsMatrix), intent(in) :: this
 end function

 !> @name Computational Methods
 !> @{

 !> @brief Returns the result of a Epetra_CrsMatrix multiplied by a Epetra_Vector x in y 
 subroutine Multiply_Vector(this,TransA,x,y,err)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   use FEpetra_Vector, only:Epetra_Vector
   class(Epetra_CrsMatrix), intent(in) :: this
   logical, intent(in) :: TransA &
   !< (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
   class(Epetra_Vector), intent(in) :: x &
   !< (In) An Epetra_Vector to multiply by.
   class(Epetra_Vector), intent(inout) :: y &
   !< (Out) An Epetra_Vector containing result. 
   type(error), optional, intent(inout) :: err &
   !< (Out) Integer error code, set to 0 if successful.
 end subroutine
 
 !> @name Computational Methods
 !> @{

 !> @brief Returns the result of a Epetra_CrsMatrix multiplied by a Epetra_MultiVector X in Y.
 subroutine Multiply_MultiVector(this,TransA,x,y,err)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   use FEpetra_MultiVector, only:Epetra_MultiVector
   class(Epetra_CrsMatrix), intent(in) :: this
   logical, intent(in) :: TransA &
   !< (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
   class(Epetra_MultiVector), intent(in) :: x &
   !< (In) An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
   class(Epetra_MultiVector), intent(inout) :: y &
   !< (Out) An Epetra_MultiVector of dimension NumVectorscontaining result.
   type(error), optional, intent(inout) :: err &
   !< (Out) Integer error code, set to 0 if successful.
 end subroutine

 !> @name Local/Global ID method
 !> @{

 !> @brief Returns true of GID is owned by the calling processor, otherwise it returns false.
 logical function MyGlobalRow(this,GID)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int), intent(in) ::GID 
 end function
 
 !> @name Attribute Accessor Methods
 !> @{

 !> @brief Returns the Epetra_Map object associated with the rows of this matrix. 
 type(Epetra_Map) function RowMatrixRowMap(this)
  use ForTrilinos_enums, only:FT_Epetra_Map_ID_t
  class(Epetra_CrsMatrix), intent(in) :: this
 end function

 !> @name Attribute Accessor Methods
 !> @{

 !> @brief Returns the Epetra_Map object associated with the rows of this matrix.
 type(Epetra_Map) function RowMap(this)
  use ForTrilinos_enums, only:FT_Epetra_Map_ID_t
  class(Epetra_CrsMatrix), intent(in) :: this
  type(FT_Epetra_Map_ID_t)            :: RowMap_ID
 end function

 !> @name Attribute Accessor Methods
 !> @{

 !> @brief Returns the current number of nonzero entries in specified global row on this processor.
 integer(c_int) function NumGlobalEntries(this,GlobalRow)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int), intent(in) :: GlobalRow
 end function

end module 

