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
! Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov) 
!*********************************************************************


#include "ForTrilinos_config.h"
#ifdef HAVE_FORTRILINOS_GALERI

module forgaleri
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  use ForTrilinos_enum_wrappers
  implicit none   ! Prevent implicit typing
#ifdef HAVE_MPI
#include <mpif.h>
#endif

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/galeri/CGaleri*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

!> @name Galeri_Utils interface
!! @{

  ! _________________ Galeri_Utils interface bodies _________________


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector* CreateCartesianCoordinates(const string CoordType, 
  !!     const Epetra_BlockMap* BlockMap, Teuchos::ParameterList& List);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Galeri_Utils_CreateCartesianCoordinates ( const char CoordType[], 
  !!     CT_Epetra_BlockMap_ID_t BlockMapID, CT_Teuchos_ParameterList_ID_t ListID );

  function Galeri_Utils_CreateCartesianCoordinates ( CoordType, BlockMapID, ListID ) result(that) &
        bind(C,name='Galeri_Utils_CreateCartesianCoordinates')
    import :: FT_Epetra_MultiVector_ID_t ,c_char ,FT_Epetra_BlockMap_ID_t , &
          FT_Teuchos_ParameterList_ID_t
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: CoordType
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: BlockMapID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ListID
  end function


  !> <BR> Original C++ prototype:
  !! void Solve(const Epetra_LinearProblem Problem);
  !> <BR> <BR> CTrilinos prototype:
  !! void Galeri_Utils_Solve_LinearProblem ( CT_Epetra_LinearProblem_ID_t ProblemID );

  subroutine Galeri_Utils_Solve_LinearProblem ( ProblemID ) &
        bind(C,name='Galeri_Utils_Solve_LinearProblem')
    import :: FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: ProblemID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void Solve(const Epetra_RowMatrix* Matrix, const Epetra_MultiVector* LHS, 
  !!     const Epetra_MultiVector* RHS);
  !> <BR> <BR> CTrilinos prototype:
  !! void Galeri_Utils_Solve_Matrix ( CT_Epetra_RowMatrix_ID_t MatrixID, 
  !!     CT_Epetra_MultiVector_ID_t LHSID, CT_Epetra_MultiVector_ID_t RHSID );

  subroutine Galeri_Utils_Solve_Matrix ( MatrixID, LHSID, RHSID ) &
        bind(C,name='Galeri_Utils_Solve_Matrix')
    import :: FT_Epetra_RowMatrix_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: MatrixID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: LHSID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: RHSID
  end subroutine


  !> <BR> Original C++ prototype:
  !! double ComputeNorm(const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS);
  !> <BR> <BR> CTrilinos prototype:
  !! double Galeri_Utils_ComputeNorm ( CT_Epetra_MultiVector_ID_t LHSID, 
  !!     CT_Epetra_MultiVector_ID_t RHSID );

  function Galeri_Utils_ComputeNorm ( LHSID, RHSID ) result(that) &
        bind(C,name='Galeri_Utils_ComputeNorm')
    import :: c_double ,FT_Epetra_MultiVector_ID_t
    
    real(c_double)                                                :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: LHSID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: RHSID
  end function


  !> <BR> Original C++ prototype:
  !! double ComputeNorm(const Epetra_RowMatrix* A, const Epetra_MultiVector* LHS, 
  !!     const Epetra_MultiVector* RHS);
  !> <BR> <BR> CTrilinos prototype:
  !! double Galeri_Utils_ComputeNorm_Matrix ( CT_Epetra_RowMatrix_ID_t AID, 
  !!     CT_Epetra_MultiVector_ID_t LHSID, CT_Epetra_MultiVector_ID_t RHSID );

  function Galeri_Utils_ComputeNorm_Matrix ( AID, LHSID, RHSID ) result(that) &
        bind(C,name='Galeri_Utils_ComputeNorm_Matrix')
    import :: c_double ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_MultiVector_ID_t
    
    real(c_double)                                                :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: LHSID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: RHSID
  end function


  !> <BR> Original C++ prototype:
  !! string toString(const int& x);
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Galeri_Utils_toString_Int ( int x );

  function Galeri_Utils_toString_Int ( x ) result(that) &
        bind(C,name='Galeri_Utils_toString_Int')
    import :: c_char ,c_int
    
    character(kind=c_char)                                        :: that
    integer(c_int)              ,intent(in)   ,value              :: x
  end function


  !> <BR> Original C++ prototype:
  !! string toString(const unsigned int& x);
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Galeri_Utils_toString_UInt ( unsigned int x );

  function Galeri_Utils_toString_UInt ( x ) result(that) &
        bind(C,name='Galeri_Utils_toString_UInt')
    import :: c_char ,c_int
    
    character(kind=c_char)                                        :: that
    integer(c_int)              ,intent(in)   ,value              :: x
  end function


  !> <BR> Original C++ prototype:
  !! string toString(const double& x);
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Galeri_Utils_toString_Double ( double x );

  function Galeri_Utils_toString_Double ( x ) result(that) &
        bind(C,name='Galeri_Utils_toString_Double')
    import :: c_char ,c_double
    
    character(kind=c_char)                                        :: that
    real(c_double)              ,intent(in)   ,value              :: x
  end function


  !> <BR> Original C++ prototype:
  !! void GetNeighboursCartesian2d(const int i, const int nx, const int ny, int & left, 
  !!     int & right, int & lower, int & upper);
  !> <BR> <BR> CTrilinos prototype:
  !! void Galeri_Utils_GetNeighboursCartesian2d ( const int i, const int nx, const int ny, 
  !!     int * left, int * right, int * lower, int * upper );

  subroutine Galeri_Utils_GetNeighboursCartesian2d ( i, nx, ny, left, right, lower, upper ) &
        bind(C,name='Galeri_Utils_GetNeighboursCartesian2d')
    import :: c_int
    
    integer(c_int)              ,intent(in)   ,value              :: i
    integer(c_int)              ,intent(in)   ,value              :: nx
    integer(c_int)              ,intent(in)   ,value              :: ny
    integer(c_int)              ,intent(inout)                    :: left
    integer(c_int)              ,intent(inout)                    :: right
    integer(c_int)              ,intent(inout)                    :: lower
    integer(c_int)              ,intent(inout)                    :: upper
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GetNeighboursCartesian2d(const int i, const int nx, const int ny, int& left, int& right, 
  !!     int& lower, int& upper, int& left2, int& right2, int& lower2, int& upper2);
  !> <BR> <BR> CTrilinos prototype:
  !! void Galeri_Utils_GetNeighboursCartesian2d_Both ( const int i, const int nx, const int ny, 
  !!     int * left, int * right, int * lower, int * upper, int * left2, int * right2, 
  !!     int * lower2, int * upper2 );

  subroutine Galeri_Utils_GetNeighboursCartesian2d_Both ( i, nx, ny, left, right, lower, &
        upper, left2, right2, lower2, upper2 ) &
        bind(C,name='Galeri_Utils_GetNeighboursCartesian2d_Both')
    import :: c_int
    
    integer(c_int)              ,intent(in)   ,value              :: i
    integer(c_int)              ,intent(in)   ,value              :: nx
    integer(c_int)              ,intent(in)   ,value              :: ny
    integer(c_int)              ,intent(inout)                    :: left
    integer(c_int)              ,intent(inout)                    :: right
    integer(c_int)              ,intent(inout)                    :: lower
    integer(c_int)              ,intent(inout)                    :: upper
    integer(c_int)              ,intent(inout)                    :: left2
    integer(c_int)              ,intent(inout)                    :: right2
    integer(c_int)              ,intent(inout)                    :: lower2
    integer(c_int)              ,intent(inout)                    :: upper2
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GetNeighboursCartesian3d(const int i, const int nx, const int ny, const int nz, 
  !!     int& left, int& right, int& lower, int& upper, int& below, int& above);
  !> <BR> <BR> CTrilinos prototype:
  !! void Galeri_Utils_GetNeighboursCartesian3d ( const int i, const int nx, const int ny, 
  !!     const int nz, int * left, int * right, int * lower, int * upper, int * below, 
  !!     int * above );

  subroutine Galeri_Utils_GetNeighboursCartesian3d ( i, nx, ny, nz, left, right, lower, &
        upper, below, above ) bind(C,name='Galeri_Utils_GetNeighboursCartesian3d')
    import :: c_int
    
    integer(c_int)              ,intent(in)   ,value              :: i
    integer(c_int)              ,intent(in)   ,value              :: nx
    integer(c_int)              ,intent(in)   ,value              :: ny
    integer(c_int)              ,intent(in)   ,value              :: nz
    integer(c_int)              ,intent(inout)                    :: left
    integer(c_int)              ,intent(inout)                    :: right
    integer(c_int)              ,intent(inout)                    :: lower
    integer(c_int)              ,intent(inout)                    :: upper
    integer(c_int)              ,intent(inout)                    :: below
    integer(c_int)              ,intent(inout)                    :: above
  end subroutine


  !> <BR> Original C++ prototype:
  !! void PrintStencil2D(const Epetra_CrsMatrix* Matrix, const int nx, const int ny, int GID = -1);
  !> <BR> <BR> CTrilinos prototype:
  !! void Galeri_Utils_PrintStencil2D ( CT_Epetra_CrsMatrix_ID_t MatrixID, const int nx, 
  !!     const int ny, int GID );

  subroutine Galeri_Utils_PrintStencil2D ( MatrixID, nx, ny, GID ) &
        bind(C,name='Galeri_Utils_PrintStencil2D')
    import :: FT_Epetra_CrsMatrix_ID_t ,c_int
    
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: MatrixID
    integer(c_int)              ,intent(in)   ,value              :: nx
    integer(c_int)              ,intent(in)   ,value              :: ny
    integer(c_int)              ,intent(in)   ,value              :: GID
  end subroutine


!> @}


!> @name Galeri_Maps interface
!! @{

  ! _________________ Galeri_Maps interface bodies _________________


  !> <BR> Original C++ prototype:
  !! Epetra_Map* CreateMap(string MapType, Epetra_Comm& Comm, Teuchos::ParameterList& List);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Galeri_Maps_CreateMap ( char MapType[], CT_Epetra_Comm_ID_t CommID, 
  !!     CT_Teuchos_ParameterList_ID_t ListID );

  function Galeri_Maps_CreateMap ( MapType, CommID, ListID ) result(that) &
        bind(C,name='Galeri_Maps_CreateMap')
    import :: FT_Epetra_Map_ID_t ,c_char ,FT_Epetra_Comm_ID_t , &
          FT_Teuchos_ParameterList_ID_t
    
    type(FT_Epetra_Map_ID_t)                                      :: that
    character(kind=c_char)                          ,dimension(*) :: MapType
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: CommID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ListID
  end function


!> @}


!> @name Galeri_CrsMatrices interface
!! @{

  ! _________________ Galeri_CrsMatrices interface bodies _________________


  !> <BR> Original C++ prototype:
  !! Epetra_CrsMatrix* CreateCrsMatrix(string MatrixType, const Epetra_Map* Map, 
  !!     Teuchos::ParameterList& List);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsMatrix_ID_t Galeri_CrsMatrices_CreateCrsMatrix ( char MatrixType[], 
  !!     CT_Epetra_Map_ID_t MapID, CT_Teuchos_ParameterList_ID_t ListID );

  function Galeri_CrsMatrices_CreateCrsMatrix ( MatrixType, MapID, ListID ) result(that) &
        bind(C,name='Galeri_CrsMatrices_CreateCrsMatrix')
    import :: FT_Epetra_CrsMatrix_ID_t ,c_char ,FT_Epetra_Map_ID_t , &
          FT_Teuchos_ParameterList_ID_t
    
    type(FT_Epetra_CrsMatrix_ID_t)                                    :: that
    character(kind=c_char)                              ,dimension(*) :: MatrixType
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: MapID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ListID
  end function


!> @}


  end interface
end module forgaleri

#endif /* HAVE_FORTRILINOS_GALERI */
