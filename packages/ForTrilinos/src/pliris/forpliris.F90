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
#ifdef HAVE_FORTRILINOS_PLIRIS

module forpliris
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  use ForTrilinos_enum_wrappers
  implicit none   ! Prevent implicit typing
#ifdef HAVE_MPI
#include <mpif.h>
#endif

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/pliris/CPliris*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

!> @name Pliris interface
!! @{

  ! _________________ Pliris interface bodies _________________





#ifdef HAVE_MPI


  !> <BR> Original C++ prototype:
  !! Pliris(Epetra_Vector * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Pliris_ID_t Pliris_Create ( CT_Epetra_Vector_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  !!     CT_Epetra_MultiVector_ID_t BID );

  function Pliris_Create ( AID, XID, BID ) result(that) bind(C,name='Pliris_Create')
    import :: FT_Pliris_ID_t ,FT_Epetra_Vector_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Pliris_ID_t)                                          :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
  end function


  !> <BR> Original C++ prototype:
  !! Pliris ();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Pliris_ID_t Pliris_Create_Default (  );

  function Pliris_Create_Default (  ) result(that) bind(C,name='Pliris_Create_Default')
    import :: FT_Pliris_ID_t
    
    type(FT_Pliris_ID_t)                                          :: that
  end function


  !> <BR> Original C++ prototype:
  !! int SetLHS(Epetra_MultiVector * X);
  !> <BR> <BR> CTrilinos prototype:
  !! int Pliris_SetLHS ( CT_Pliris_ID_t selfID, CT_Epetra_MultiVector_ID_t XID );

  function Pliris_SetLHS ( selfID, XID ) result(that) bind(C,name='Pliris_SetLHS')
    import :: c_int ,FT_Pliris_ID_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Pliris_ID_t)        ,intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
  end function


  !> <BR> Original C++ prototype:
  !! int SetRHS(Epetra_MultiVector * B);
  !> <BR> <BR> CTrilinos prototype:
  !! int Pliris_SetRHS ( CT_Pliris_ID_t selfID, CT_Epetra_MultiVector_ID_t BID );

  function Pliris_SetRHS ( selfID, BID ) result(that) bind(C,name='Pliris_SetRHS')
    import :: c_int ,FT_Pliris_ID_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Pliris_ID_t)        ,intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
  end function


  !> <BR> Original C++ prototype:
  !! int SetMatrix(Epetra_Vector * A);
  !> <BR> <BR> CTrilinos prototype:
  !! int Pliris_SetMatrix ( CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID );

  function Pliris_SetMatrix ( selfID, AID ) result(that) bind(C,name='Pliris_SetMatrix')
    import :: c_int ,FT_Pliris_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Pliris_ID_t)        ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: AID
  end function


  !> <BR> Original C++ prototype:
  !! int SetMatrix(Epetra_SerialDenseVector * A);
  !> <BR> <BR> CTrilinos prototype:
  !! int Pliris_SetMatrix_Serial ( CT_Pliris_ID_t selfID, CT_Epetra_SerialDenseVector_ID_t AID );

  function Pliris_SetMatrix_Serial ( selfID, AID ) result(that) &
        bind(C,name='Pliris_SetMatrix_Serial')
    import :: c_int ,FT_Pliris_ID_t ,FT_Epetra_SerialDenseVector_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Pliris_ID_t)        ,intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: AID
  end function


  !> <BR> Original C++ prototype:
  !! int GetDistribution( int * nprocs_row, int * number_of_unknowns, int * nrhs, int * my_rows, 
  !!     int * my_cols, int * my_first_row, int * my_first_col, int * my_rhs, int * my_row, 
  !!     int * my_col );
  !> <BR> <BR> CTrilinos prototype:
  !! int Pliris_GetDistribution ( CT_Pliris_ID_t selfID, int * nprocs_row, int * number_of_unknowns, 
  !!     int * nrhs, int * my_rows, int * my_cols, int * my_first_row, int * my_first_col, 
  !!     int * my_rhs, int * my_row, int * my_col );

  function Pliris_GetDistribution ( selfID, nprocs_row, number_of_unknowns, nrhs, my_rows, &
        my_cols, my_first_row, my_first_col, my_rhs, my_row, my_col ) result(that) &
        bind(C,name='Pliris_GetDistribution')
    import :: c_int ,FT_Pliris_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Pliris_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)       ,intent(in)                       :: nprocs_row
    integer(c_int)       ,intent(in)                       :: number_of_unknowns
    integer(c_int)       ,intent(in)                       :: nrhs
    integer(c_int)       ,intent(out)                      :: my_rows
    integer(c_int)       ,intent(out)                      :: my_cols
    integer(c_int)       ,intent(out)                      :: my_first_row
    integer(c_int)       ,intent(out)                      :: my_first_col
    integer(c_int)       ,intent(out)                      :: my_rhs
    integer(c_int)       ,intent(out)                      :: my_row
    integer(c_int)       ,intent(out)                      :: my_col
  end function


  !> <BR> Original C++ prototype:
  !! int FactorSolve( Epetra_Vector * A, int my_rows, int my_cols, int* matrix_size, int* num_procsr, 
  !!     int* num_rhs, double* secs);
  !> <BR> <BR> CTrilinos prototype:
  !! int Pliris_FactorSolve ( CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID, int my_rows, 
  !!     int my_cols, int * matrix_size, int * num_procsr, int * num_rhs, double * secs );

  function Pliris_FactorSolve ( selfID, AID, my_rows, my_cols, matrix_size, num_procsr, &
        num_rhs, secs ) result(that) bind(C,name='Pliris_FactorSolve')
    import :: c_int ,FT_Pliris_ID_t ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Pliris_ID_t)        ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: AID
    integer(c_int)              ,intent(in)   ,value              :: my_rows
    integer(c_int)              ,intent(in)   ,value              :: my_cols
    integer(c_int)              ,intent(in)                       :: matrix_size
    integer(c_int)              ,intent(in)                       :: num_procsr
    integer(c_int)              ,intent(in)                       :: num_rhs
    real(c_double)              ,intent(out)                      :: secs
  end function


  !> <BR> Original C++ prototype:
  !! int FactorSolve( Epetra_SerialDenseVector * AA, int my_rows, int my_cols, int* matrix_size, 
  !!     int* num_procsr, int* num_rhs, double* secs);
  !> <BR> <BR> CTrilinos prototype:
  !! int Pliris_FactorSolve_Serial ( CT_Pliris_ID_t selfID, CT_Epetra_SerialDenseVector_ID_t AAID, 
  !!     int my_rows, int my_cols, int * matrix_size, int * num_procsr, int * num_rhs, 
  !!     double * secs );

! function Pliris_FactorSolve_Serial ( selfID, AAID, my_rows, my_cols, matrix_size, &
!       num_procsr, num_rhs, secs ) result(that) bind(C,name='Pliris_FactorSolve_Serial')
!   import :: c_int ,FT_Pliris_ID_t ,FT_Epetra_SerialDenseVector_ID_t ,c_double
!   
!   integer(c_int)                                                :: that
!   type(FT_Pliris_ID_t)        ,intent(in)   ,value              :: selfID
!   type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: AAID
!   integer(c_int)              ,intent(in)   ,value              :: my_rows
!   integer(c_int)              ,intent(in)   ,value              :: my_cols
!   integer(c_int)              ,intent(in)   ,value              :: matrix_size
!   integer(c_int)              ,intent(in)   ,value              :: num_procsr
!   integer(c_int)              ,intent(in)   ,value              :: num_rhs
!   real(c_double)              ,intent(out)                      :: secs
! end function


  !> <BR> Original C++ prototype:
  !! int Factor( Epetra_Vector* A, int* matrix_size, int* num_procsr, int* permute, double* secs);
  !> <BR> <BR> CTrilinos prototype:
  !! int Pliris_Factor ( CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID, int * matrix_size, 
  !!     int * num_procsr, int * permute, double * secs );

  function Pliris_Factor ( selfID, AID, matrix_size, num_procsr, permute, secs ) &
        result(that) bind(C,name='Pliris_Factor')
    import :: c_int ,FT_Pliris_ID_t ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Pliris_ID_t)        ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: AID
    integer(c_int)              ,intent(in)   ,value              :: matrix_size
    integer(c_int)              ,intent(in)   ,value              :: num_procsr
    integer(c_int)              ,intent(out)        ,dimension(*) :: permute
    real(c_double)              ,intent(out)                      :: secs
  end function


  !> <BR> Original C++ prototype:
  !! int Solve(int* permute, int* num_rhs);
  !> <BR> <BR> CTrilinos prototype:
  !! int Pliris_Solve ( CT_Pliris_ID_t selfID, int * permute, int * num_rhs );

  function Pliris_Solve ( selfID, permute, num_rhs ) result(that) &
        bind(C,name='Pliris_Solve')
    import :: c_int ,FT_Pliris_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Pliris_ID_t)        ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)         ,dimension(*) :: permute
    integer(c_int)              ,intent(in)   ,value              :: num_rhs
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Pliris(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Pliris_Destroy ( CT_Pliris_ID_t * selfID );

  subroutine Pliris_Destroy ( selfID ) bind(C,name='Pliris_Destroy')
    import :: FT_Pliris_ID_t
    
    type(FT_Pliris_ID_t)                                          :: selfID
  end subroutine


#endif /* HAVE_MPI */


!> @}


  end interface
end module forpliris

#endif
