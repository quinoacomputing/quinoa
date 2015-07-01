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
#ifdef HAVE_FORTRILINOS_AZTECOO

module foraztecoo
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  use ForTrilinos_enum_wrappers
  implicit none   ! Prevent implicit typing
#ifdef HAVE_MPI
#include <mpif.h>
#endif

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/aztecoo/CAztecoo*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

!> @name AztecOO interface
!! @{

  ! _________________ AztecOO interface bodies _________________


  !> <BR> Original C++ prototype:
  !! AztecOO(Epetra_Operator * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_ID_t AztecOO_Create_FromOperator ( CT_Epetra_Operator_ID_t AID, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID );

  function AztecOO_Create_FromOperator ( AID, XID, BID ) result(that) &
        bind(C,name='AztecOO_Create_FromOperator')
    import :: FT_AztecOO_ID_t ,FT_Epetra_Operator_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_AztecOO_ID_t)                                         :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_ID_t AztecOO_Create_FromRowMatrix ( CT_Epetra_RowMatrix_ID_t AID, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID );

  function AztecOO_Create_FromRowMatrix ( AID, XID, BID ) result(that) &
        bind(C,name='AztecOO_Create_FromRowMatrix')
    import :: FT_AztecOO_ID_t ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_AztecOO_ID_t)                                         :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO(const Epetra_LinearProblem& LinearProblem);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_ID_t AztecOO_Create_FromLinearProblem ( CT_Epetra_LinearProblem_ID_t LinearProblemID );

  function AztecOO_Create_FromLinearProblem ( LinearProblemID ) result(that) &
        bind(C,name='AztecOO_Create_FromLinearProblem')
    import :: FT_AztecOO_ID_t ,FT_Epetra_LinearProblem_ID_t
    
    type(FT_AztecOO_ID_t)                                         :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: LinearProblemID
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_ID_t AztecOO_Create (  );

  function AztecOO_Create (  ) result(that) bind(C,name='AztecOO_Create')
    import :: FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)                                         :: that
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO(const AztecOO& Solver);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_ID_t AztecOO_Duplicate ( CT_AztecOO_ID_t SolverID );

  function AztecOO_Duplicate ( SolverID ) result(that) bind(C,name='AztecOO_Duplicate')
    import :: FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)                                         :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: SolverID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~AztecOO(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void AztecOO_Destroy ( CT_AztecOO_ID_t * selfID );

  subroutine AztecOO_Destroy ( selfID ) bind(C,name='AztecOO_Destroy')
    import :: FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)                                         :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int SetProblem(const Epetra_LinearProblem& prob, bool call_SetPrecMatrix=false);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetProblem ( CT_AztecOO_ID_t selfID, CT_Epetra_LinearProblem_ID_t probID, 
  !!     boolean call_SetPrecMatrix );

  function AztecOO_SetProblem ( selfID, probID, call_SetPrecMatrix ) result(that) &
        bind(C,name='AztecOO_SetProblem')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_LinearProblem_ID_t ,FT_boolean_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: probID
    integer(FT_boolean_t)       ,intent(in)   ,value              :: call_SetPrecMatrix
  end function


  !> <BR> Original C++ prototype:
  !! int SetUserOperator(Epetra_Operator * UserOperator);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetUserOperator ( CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t UserOperatorID );

  function AztecOO_SetUserOperator ( selfID, UserOperatorID ) result(that) &
        bind(C,name='AztecOO_SetUserOperator')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_Operator_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: UserOperatorID
  end function


  !> <BR> Original C++ prototype:
  !! int SetUserMatrix(Epetra_RowMatrix * UserMatrix, bool call_SetPrecMatrix=false);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetUserMatrix ( CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t UserMatrixID, 
  !!     boolean call_SetPrecMatrix );

  function AztecOO_SetUserMatrix ( selfID, UserMatrixID, call_SetPrecMatrix ) result(that) &
        bind(C,name='AztecOO_SetUserMatrix')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_RowMatrix_ID_t ,FT_boolean_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: UserMatrixID
    integer(FT_boolean_t)       ,intent(in)   ,value              :: call_SetPrecMatrix
  end function


  !> <BR> Original C++ prototype:
  !! int SetLHS(Epetra_MultiVector * X);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetLHS ( CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t XID );

  function AztecOO_SetLHS ( selfID, XID ) result(that) bind(C,name='AztecOO_SetLHS')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
  end function


  !> <BR> Original C++ prototype:
  !! int SetRHS(Epetra_MultiVector * B);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetRHS ( CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t BID );

  function AztecOO_SetRHS ( selfID, BID ) result(that) bind(C,name='AztecOO_SetRHS')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
  end function


  !> <BR> Original C++ prototype:
  !! int SetPrecMatrix(Epetra_RowMatrix * PrecMatrix);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetPrecMatrix ( CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t PrecMatrixID );

  function AztecOO_SetPrecMatrix ( selfID, PrecMatrixID ) result(that) &
        bind(C,name='AztecOO_SetPrecMatrix')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: PrecMatrixID
  end function


  !> <BR> Original C++ prototype:
  !! int SetPrecOperator(Epetra_Operator * PrecOperator);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetPrecOperator ( CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t PrecOperatorID );

  function AztecOO_SetPrecOperator ( selfID, PrecOperatorID ) result(that) &
        bind(C,name='AztecOO_SetPrecOperator')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_Operator_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: PrecOperatorID
  end function


  !> <BR> Original C++ prototype:
  !! int ConstructPreconditioner(double & condest);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_ConstructPreconditioner ( CT_AztecOO_ID_t selfID, double * condest );

  function AztecOO_ConstructPreconditioner ( selfID, condest ) result(that) &
        bind(C,name='AztecOO_ConstructPreconditioner')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    real(c_double)              ,intent(inout)                    :: condest
  end function


  !> <BR> Original C++ prototype:
  !! int DestroyPreconditioner();
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_DestroyPreconditioner ( CT_AztecOO_ID_t selfID );

  function AztecOO_DestroyPreconditioner ( selfID ) result(that) &
        bind(C,name='AztecOO_DestroyPreconditioner')
    import :: c_int ,FT_AztecOO_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double Condest() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double AztecOO_Condest ( CT_AztecOO_ID_t selfID );

  function AztecOO_Condest ( selfID ) result(that) bind(C,name='AztecOO_Condest')
    import :: c_double ,FT_AztecOO_ID_t
    
    real(c_double)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int CheckInput() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_CheckInput ( CT_AztecOO_ID_t selfID );

  function AztecOO_CheckInput ( selfID ) result(that) bind(C,name='AztecOO_CheckInput')
    import :: c_int ,FT_AztecOO_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_LinearProblem * GetProblem() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_LinearProblem_ID_t AztecOO_GetProblem ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetProblem ( selfID ) result(that) bind(C,name='AztecOO_GetProblem')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_AztecOO_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t)                                  :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Operator * GetUserOperator() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Operator_ID_t AztecOO_GetUserOperator ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetUserOperator ( selfID ) result(that) &
        bind(C,name='AztecOO_GetUserOperator')
    import :: FT_Epetra_Operator_ID_t ,FT_AztecOO_ID_t
    
    type(FT_Epetra_Operator_ID_t)                                  :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_RowMatrix * GetUserMatrix() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_RowMatrix_ID_t AztecOO_GetUserMatrix ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetUserMatrix ( selfID ) result(that) &
        bind(C,name='AztecOO_GetUserMatrix')
    import :: FT_Epetra_RowMatrix_ID_t ,FT_AztecOO_ID_t
    
    type(FT_Epetra_RowMatrix_ID_t)                                  :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Operator * GetPrecOperator() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Operator_ID_t AztecOO_GetPrecOperator ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetPrecOperator ( selfID ) result(that) &
        bind(C,name='AztecOO_GetPrecOperator')
    import :: FT_Epetra_Operator_ID_t ,FT_AztecOO_ID_t
    
    type(FT_Epetra_Operator_ID_t)                                  :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_RowMatrix * GetPrecMatrix() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_RowMatrix_ID_t AztecOO_GetPrecMatrix ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetPrecMatrix ( selfID ) result(that) &
        bind(C,name='AztecOO_GetPrecMatrix')
    import :: FT_Epetra_RowMatrix_ID_t ,FT_AztecOO_ID_t
    
    type(FT_Epetra_RowMatrix_ID_t)                                  :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector * GetLHS() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t AztecOO_GetLHS ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetLHS ( selfID ) result(that) bind(C,name='AztecOO_GetLHS')
    import :: FT_Epetra_MultiVector_ID_t ,FT_AztecOO_ID_t
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector * GetRHS() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t AztecOO_GetRHS ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetRHS ( selfID ) result(that) bind(C,name='AztecOO_GetRHS')
    import :: FT_Epetra_MultiVector_ID_t ,FT_AztecOO_ID_t
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void PrintLinearSystem(const char* name);
  !> <BR> <BR> CTrilinos prototype:
  !! void AztecOO_PrintLinearSystem ( CT_AztecOO_ID_t selfID, const char * name );

  subroutine AztecOO_PrintLinearSystem ( selfID, name ) &
        bind(C,name='AztecOO_PrintLinearSystem')
    import :: FT_AztecOO_ID_t ,c_char
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: name
  end subroutine


  !> <BR> Original C++ prototype:
  !! int SetParameters(Teuchos::ParameterList& parameterlist, bool cerr_warning_if_unused=false);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetParameters ( CT_AztecOO_ID_t selfID, 
  !!     CT_Teuchos_ParameterList_ID_t parameterlistID, boolean cerr_warning_if_unused );

  function AztecOO_SetParameters ( selfID, parameterlistID, cerr_warning_if_unused ) result(that) &
        bind(C,name='AztecOO_SetParameters')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Teuchos_ParameterList_ID_t ,FT_boolean_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: parameterlistID
    integer(FT_boolean_t)       ,intent(in)   ,value              :: cerr_warning_if_unused
  end function


  !> <BR> Original C++ prototype:
  !! int SetAztecDefaults();
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetAztecDefaults ( CT_AztecOO_ID_t selfID );

  function AztecOO_SetAztecDefaults ( selfID ) result(that) &
        bind(C,name='AztecOO_SetAztecDefaults')
    import :: c_int ,FT_AztecOO_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int SetAztecOption(int option, int value);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetAztecOption ( CT_AztecOO_ID_t selfID, int option, int value );

  function AztecOO_SetAztecOption ( selfID, option, value ) result(that) &
        bind(C,name='AztecOO_SetAztecOption')
    import :: c_int ,FT_AztecOO_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: option
    integer(c_int)              ,intent(in)   ,value              :: value
  end function


  !> <BR> Original C++ prototype:
  !! int GetAztecOption(int option);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_GetAztecOption ( CT_AztecOO_ID_t selfID, int option );

  function AztecOO_GetAztecOption ( selfID, option ) result(that) &
        bind(C,name='AztecOO_GetAztecOption')
    import :: c_int ,FT_AztecOO_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: option
  end function


  !> <BR> Original C++ prototype:
  !! int SetAztecParam(int param, double value);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetAztecParam ( CT_AztecOO_ID_t selfID, int param, double value );

  function AztecOO_SetAztecParam ( selfID, param, value ) result(that) &
        bind(C,name='AztecOO_SetAztecParam')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: param
    real(c_double)              ,intent(in)   ,value              :: value
  end function


  !> <BR> Original C++ prototype:
  !! const int* GetAllAztecOptions() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const int * AztecOO_GetAllAztecOptions ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetAllAztecOptions ( selfID ) result(that) &
        bind(C,name='AztecOO_GetAllAztecOptions')
    import :: c_ptr ,FT_AztecOO_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const double* GetAllAztecParams() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const double * AztecOO_GetAllAztecParams ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetAllAztecParams ( selfID ) result(that) &
        bind(C,name='AztecOO_GetAllAztecParams')
    import :: c_ptr ,FT_AztecOO_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int SetAllAztecOptions(const int * options);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetAllAztecOptions ( CT_AztecOO_ID_t selfID, const int * options );

  function AztecOO_SetAllAztecOptions ( selfID, options ) result(that) &
        bind(C,name='AztecOO_SetAllAztecOptions')
    import :: c_int ,FT_AztecOO_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)         ,dimension(*) :: options
  end function


  !> <BR> Original C++ prototype:
  !! int SetAllAztecParams(const double * params);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetAllAztecParams ( CT_AztecOO_ID_t selfID, const double * params );

  function AztecOO_SetAllAztecParams ( selfID, params ) result(that) &
        bind(C,name='AztecOO_SetAllAztecParams')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    real(c_double)              ,intent(in)         ,dimension(*) :: params
  end function


  !> <BR> Original C++ prototype:
  !! int Iterate(int MaxIters, double Tolerance);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_Iterate_Current ( CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance );

  function AztecOO_Iterate_Current ( selfID, MaxIters, Tolerance ) result(that) &
        bind(C,name='AztecOO_Iterate_Current')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: MaxIters
    real(c_double)              ,intent(in)   ,value              :: Tolerance
  end function


  !> <BR> Original C++ prototype:
  !! int Iterate(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B, 
  !!     int MaxIters, double Tolerance);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_Iterate ( CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID, int MaxIters, 
  !!     double Tolerance );

  function AztecOO_Iterate ( selfID, AID, XID, BID, MaxIters, Tolerance ) result(that) &
        bind(C,name='AztecOO_Iterate')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_MultiVector_ID_t , &
          c_double
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
    integer(c_int)              ,intent(in)   ,value              :: MaxIters
    real(c_double)              ,intent(in)   ,value              :: Tolerance
  end function


  !> <BR> Original C++ prototype:
  !! int recursiveIterate(int MaxIters, double Tolerance);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_recursiveIterate ( CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance );

  function AztecOO_recursiveIterate ( selfID, MaxIters, Tolerance ) result(that) &
        bind(C,name='AztecOO_recursiveIterate')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: MaxIters
    real(c_double)              ,intent(in)   ,value              :: Tolerance
  end function


  !> <BR> Original C++ prototype:
  !! const double *GetAztecStatus() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const double * AztecOO_GetAztecStatus ( CT_AztecOO_ID_t selfID );

  function AztecOO_GetAztecStatus ( selfID ) result(that) &
        bind(C,name='AztecOO_GetAztecStatus')
    import :: c_ptr ,FT_AztecOO_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int SetUseAdaptiveDefaultsTrue();
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetUseAdaptiveDefaultsTrue ( CT_AztecOO_ID_t selfID );

  function AztecOO_SetUseAdaptiveDefaultsTrue ( selfID ) result(that) &
        bind(C,name='AztecOO_SetUseAdaptiveDefaultsTrue')
    import :: c_int ,FT_AztecOO_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int SetAdaptiveParams(int NumTrials, double * athresholds, double * rthresholds, 
  !!     double condestThreshold, double maxFill, int maxKspace);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_SetAdaptiveParams ( CT_AztecOO_ID_t selfID, int NumTrials, double * athresholds, 
  !!     double * rthresholds, double condestThreshold, double maxFill, int maxKspace );

  function AztecOO_SetAdaptiveParams ( selfID, NumTrials, athresholds, rthresholds, &
        condestThreshold, maxFill, maxKspace ) result(that) &
        bind(C,name='AztecOO_SetAdaptiveParams')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumTrials
    real(c_double)                                  ,dimension(*) :: athresholds
    real(c_double)                                  ,dimension(*) :: rthresholds
    real(c_double)              ,intent(in)   ,value              :: condestThreshold
    real(c_double)              ,intent(in)   ,value              :: maxFill
    integer(c_int)              ,intent(in)   ,value              :: maxKspace
  end function


  !> <BR> Original C++ prototype:
  !! int AdaptiveIterate(int MaxIters, int MaxSolveAttempts, double Tolerance);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_AdaptiveIterate ( CT_AztecOO_ID_t selfID, int MaxIters, int MaxSolveAttempts, 
  !!     double Tolerance );

  function AztecOO_AdaptiveIterate ( selfID, MaxIters, MaxSolveAttempts, Tolerance ) result(that) &
        bind(C,name='AztecOO_AdaptiveIterate')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: MaxIters
    integer(c_int)              ,intent(in)   ,value              :: MaxSolveAttempts
    real(c_double)              ,intent(in)   ,value              :: Tolerance
  end function


  !> <BR> Original C++ prototype:
  !! int NumIters() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_NumIters ( CT_AztecOO_ID_t selfID );

  function AztecOO_NumIters ( selfID ) result(that) bind(C,name='AztecOO_NumIters')
    import :: c_int ,FT_AztecOO_ID_t
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double TrueResidual() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double AztecOO_TrueResidual ( CT_AztecOO_ID_t selfID );

  function AztecOO_TrueResidual ( selfID ) result(that) bind(C,name='AztecOO_TrueResidual')
    import :: c_double ,FT_AztecOO_ID_t
    
    real(c_double)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double ScaledResidual() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double AztecOO_ScaledResidual ( CT_AztecOO_ID_t selfID );

  function AztecOO_ScaledResidual ( selfID ) result(that) &
        bind(C,name='AztecOO_ScaledResidual')
    import :: c_double ,FT_AztecOO_ID_t
    
    real(c_double)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double RecursiveResidual() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double AztecOO_RecursiveResidual ( CT_AztecOO_ID_t selfID );

  function AztecOO_RecursiveResidual ( selfID ) result(that) &
        bind(C,name='AztecOO_RecursiveResidual')
    import :: c_double ,FT_AztecOO_ID_t
    
    real(c_double)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double SolveTime() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double AztecOO_SolveTime ( CT_AztecOO_ID_t selfID );

  function AztecOO_SolveTime ( selfID ) result(that) bind(C,name='AztecOO_SolveTime')
    import :: c_double ,FT_AztecOO_ID_t
    
    real(c_double)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int GetAllAztecStatus(double * status);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_GetAllAztecStatus ( CT_AztecOO_ID_t selfID, double * status );

  function AztecOO_GetAllAztecStatus ( selfID, status ) result(that) &
        bind(C,name='AztecOO_GetAllAztecStatus')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: status
  end function


!> @}


!> @name AztecOO_StatusTest interface
!! @{

  ! _________________ AztecOO_StatusTest interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest_Degeneralize ( CTrilinos_Universal_ID_t id );

  function AztecOO_StatusTest_Degeneralize ( id ) result(that) &
        bind(C,name='AztecOO_StatusTest_Degeneralize')
    import :: FT_AztecOO_StatusTest_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_AztecOO_StatusTest_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t AztecOO_StatusTest_Generalize ( CT_AztecOO_StatusTest_ID_t id );

  function AztecOO_StatusTest_Generalize ( id ) result(that) &
        bind(C,name='AztecOO_StatusTest_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_AztecOO_StatusTest_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_AztecOO_StatusTest_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~AztecOO_StatusTest();
  !> <BR> <BR> CTrilinos prototype:
  !! void AztecOO_StatusTest_Destroy ( CT_AztecOO_StatusTest_ID_t * selfID );

  subroutine AztecOO_StatusTest_Destroy ( selfID ) &
        bind(C,name='AztecOO_StatusTest_Destroy')
    import :: FT_AztecOO_StatusTest_ID_t
    
    type(FT_AztecOO_StatusTest_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual bool ResidualVectorRequired() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean AztecOO_StatusTest_ResidualVectorRequired ( CT_AztecOO_StatusTest_ID_t selfID );

  function AztecOO_StatusTest_ResidualVectorRequired ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTest_ResidualVectorRequired')
    import :: FT_boolean_t ,FT_AztecOO_StatusTest_ID_t
    
    integer(FT_boolean_t)                                             :: that
    type(FT_AztecOO_StatusTest_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
  !!     double CurrentResNormEst, bool SolutionUpdated) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusType_E_t AztecOO_StatusTest_CheckStatus ( CT_AztecOO_StatusTest_ID_t selfID, 
  !!     int CurrentIter, CT_Epetra_MultiVector_ID_t CurrentResVectorID, double CurrentResNormEst, 
  !!     boolean SolutionUpdated );

  function AztecOO_StatusTest_CheckStatus ( selfID, CurrentIter, CurrentResVectorID, &
        CurrentResNormEst, SolutionUpdated ) result(that) &
        bind(C,name='AztecOO_StatusTest_CheckStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTest_ID_t ,c_int , &
          FT_Epetra_MultiVector_ID_t ,c_double ,FT_boolean_t
    
    integer(FT_AztecOO_StatusType_E_t)                                  :: that
    type(FT_AztecOO_StatusTest_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: CurrentIter
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: CurrentResVectorID
    real(c_double)                  ,intent(in)   ,value              :: CurrentResNormEst
    integer(FT_boolean_t)           ,intent(in)   ,value              :: SolutionUpdated
  end function


  !> <BR> Original C++ prototype:
  !! virtual AztecOO_StatusType GetStatus() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusType_E_t AztecOO_StatusTest_GetStatus ( CT_AztecOO_StatusTest_ID_t selfID );

  function AztecOO_StatusTest_GetStatus ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTest_GetStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTest_ID_t
    
    integer(FT_AztecOO_StatusType_E_t)                                  :: that
    type(FT_AztecOO_StatusTest_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name AztecOO_StatusTestCombo interface
!! @{

  ! _________________ AztecOO_StatusTestCombo interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Degeneralize ( CTrilinos_Universal_ID_t id );

  function AztecOO_StatusTestCombo_Degeneralize ( id ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_Degeneralize')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t)     ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t AztecOO_StatusTestCombo_Generalize ( CT_AztecOO_StatusTestCombo_ID_t id );

  function AztecOO_StatusTestCombo_Generalize ( id ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_AztecOO_StatusTestCombo_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                       :: that
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusTestCombo(ComboType t);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create ( CT_ComboType_E_t t );

  function AztecOO_StatusTestCombo_Create ( t ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_Create')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,FT_ComboType_E_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t)                                  :: that
    integer(FT_ComboType_E_t)            ,intent(in)   ,value              :: t
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_OneTest ( CT_ComboType_E_t t, 
  !!     CT_AztecOO_StatusTest_ID_t aID );

  function AztecOO_StatusTestCombo_Create_OneTest ( t, aID ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_Create_OneTest')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,FT_ComboType_E_t , &
          FT_AztecOO_StatusTest_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t)                                  :: that
    integer(FT_ComboType_E_t)            ,intent(in)   ,value              :: t
    type(FT_AztecOO_StatusTest_ID_t)     ,intent(in)   ,value              :: aID
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a, AztecOO_StatusTest& b);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_TwoTests ( CT_ComboType_E_t t, 
  !!     CT_AztecOO_StatusTest_ID_t aID, CT_AztecOO_StatusTest_ID_t bID );

  function AztecOO_StatusTestCombo_Create_TwoTests ( t, aID, bID ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_Create_TwoTests')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,FT_ComboType_E_t , &
          FT_AztecOO_StatusTest_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t)                                  :: that
    integer(FT_ComboType_E_t)            ,intent(in)   ,value              :: t
    type(FT_AztecOO_StatusTest_ID_t)     ,intent(in)   ,value              :: aID
    type(FT_AztecOO_StatusTest_ID_t)     ,intent(in)   ,value              :: bID
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusTestCombo& AddStatusTest(AztecOO_StatusTest& a);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_AddStatusTest ( CT_AztecOO_StatusTestCombo_ID_t selfID, 
  !!     CT_AztecOO_StatusTest_ID_t aID );

  function AztecOO_StatusTestCombo_AddStatusTest ( selfID, aID ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_AddStatusTest')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,FT_AztecOO_StatusTest_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t)                                  :: that
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
    type(FT_AztecOO_StatusTest_ID_t)     ,intent(in)   ,value              :: aID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~AztecOO_StatusTestCombo();
  !> <BR> <BR> CTrilinos prototype:
  !! void AztecOO_StatusTestCombo_Destroy ( CT_AztecOO_StatusTestCombo_ID_t * selfID );

  subroutine AztecOO_StatusTestCombo_Destroy ( selfID ) &
        bind(C,name='AztecOO_StatusTestCombo_Destroy')
    import :: FT_AztecOO_StatusTestCombo_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! bool ResidualVectorRequired() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean AztecOO_StatusTestCombo_ResidualVectorRequired ( CT_AztecOO_StatusTestCombo_ID_t selfID );

  function AztecOO_StatusTestCombo_ResidualVectorRequired ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_ResidualVectorRequired')
    import :: FT_boolean_t ,FT_AztecOO_StatusTestCombo_ID_t
    
    integer(FT_boolean_t)                                                  :: that
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
  !!     double CurrentResNormEst, bool SolutionUpdated);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_CheckStatus ( CT_AztecOO_StatusTestCombo_ID_t selfID, 
  !!     int CurrentIter, CT_Epetra_MultiVector_ID_t CurrentResVectorID, double CurrentResNormEst, 
  !!     boolean SolutionUpdated );

  function AztecOO_StatusTestCombo_CheckStatus ( selfID, CurrentIter, CurrentResVectorID, &
        CurrentResNormEst, SolutionUpdated ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_CheckStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestCombo_ID_t ,c_int , &
          FT_Epetra_MultiVector_ID_t ,c_double ,FT_boolean_t
    
    integer(FT_AztecOO_StatusType_E_t)                                     :: that
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                       ,intent(in)   ,value              :: CurrentIter
    type(FT_Epetra_MultiVector_ID_t)     ,intent(in)   ,value              :: CurrentResVectorID
    real(c_double)                       ,intent(in)   ,value              :: CurrentResNormEst
    integer(FT_boolean_t)                ,intent(in)   ,value              :: SolutionUpdated
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusType GetStatus() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_GetStatus ( CT_AztecOO_StatusTestCombo_ID_t selfID );

  function AztecOO_StatusTestCombo_GetStatus ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_GetStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestCombo_ID_t
    
    integer(FT_AztecOO_StatusType_E_t)                                     :: that
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! ComboType GetComboType() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_ComboType_E_t AztecOO_StatusTestCombo_GetComboType ( CT_AztecOO_StatusTestCombo_ID_t selfID );

  function AztecOO_StatusTestCombo_GetComboType ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestCombo_GetComboType')
    import :: FT_ComboType_E_t ,FT_AztecOO_StatusTestCombo_ID_t
    
    integer(FT_ComboType_E_t)                                              :: that
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name AztecOO_StatusTestMaxIters interface
!! @{

  ! _________________ AztecOO_StatusTestMaxIters interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Degeneralize ( CTrilinos_Universal_ID_t id );

  function AztecOO_StatusTestMaxIters_Degeneralize ( id ) result(that) &
        bind(C,name='AztecOO_StatusTestMaxIters_Degeneralize')
    import :: FT_AztecOO_StatusTestMaxIters_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t)        ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t AztecOO_StatusTestMaxIters_Generalize ( CT_AztecOO_StatusTestMaxIters_ID_t id );

  function AztecOO_StatusTestMaxIters_Generalize ( id ) result(that) &
        bind(C,name='AztecOO_StatusTestMaxIters_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                          :: that
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusTestMaxIters(int MaxIters);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Create ( int MaxIters );

  function AztecOO_StatusTestMaxIters_Create ( MaxIters ) result(that) &
        bind(C,name='AztecOO_StatusTestMaxIters_Create')
    import :: FT_AztecOO_StatusTestMaxIters_ID_t ,c_int
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t)                                  :: that
    integer(c_int)                          ,intent(in)   ,value              :: MaxIters
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~AztecOO_StatusTestMaxIters();
  !> <BR> <BR> CTrilinos prototype:
  !! void AztecOO_StatusTestMaxIters_Destroy ( CT_AztecOO_StatusTestMaxIters_ID_t * selfID );

  subroutine AztecOO_StatusTestMaxIters_Destroy ( selfID ) &
        bind(C,name='AztecOO_StatusTestMaxIters_Destroy')
    import :: FT_AztecOO_StatusTestMaxIters_ID_t
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! bool ResidualVectorRequired() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean AztecOO_StatusTestMaxIters_ResidualVectorRequired ( CT_AztecOO_StatusTestMaxIters_ID_t selfID );

  function AztecOO_StatusTestMaxIters_ResidualVectorRequired ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestMaxIters_ResidualVectorRequired')
    import :: FT_boolean_t ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    integer(FT_boolean_t)                                                     :: that
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
  !!     double CurrentResNormEst, bool SolutionUpdated);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusType_E_t AztecOO_StatusTestMaxIters_CheckStatus ( CT_AztecOO_StatusTestMaxIters_ID_t selfID, 
  !!     int CurrentIter, CT_Epetra_MultiVector_ID_t CurrentResVectorID, double CurrentResNormEst, 
  !!     boolean SolutionUpdated );

  function AztecOO_StatusTestMaxIters_CheckStatus ( selfID, CurrentIter, CurrentResVectorID, &
        CurrentResNormEst, SolutionUpdated ) result(that) &
        bind(C,name='AztecOO_StatusTestMaxIters_CheckStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestMaxIters_ID_t ,c_int , &
          FT_Epetra_MultiVector_ID_t ,c_double ,FT_boolean_t
    
    integer(FT_AztecOO_StatusType_E_t)                                        :: that
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                          ,intent(in)   ,value              :: CurrentIter
    type(FT_Epetra_MultiVector_ID_t)        ,intent(in)   ,value              :: CurrentResVectorID
    real(c_double)                          ,intent(in)   ,value              :: CurrentResNormEst
    integer(FT_boolean_t)                   ,intent(in)   ,value              :: SolutionUpdated
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusType GetStatus() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusType_E_t AztecOO_StatusTestMaxIters_GetStatus ( CT_AztecOO_StatusTestMaxIters_ID_t selfID );

  function AztecOO_StatusTestMaxIters_GetStatus ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestMaxIters_GetStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    integer(FT_AztecOO_StatusType_E_t)                                        :: that
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int GetMaxIters() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_StatusTestMaxIters_GetMaxIters ( CT_AztecOO_StatusTestMaxIters_ID_t selfID );

  function AztecOO_StatusTestMaxIters_GetMaxIters ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestMaxIters_GetMaxIters')
    import :: c_int ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    integer(c_int)                                                            :: that
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int GetNumIters() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_StatusTestMaxIters_GetNumIters ( CT_AztecOO_StatusTestMaxIters_ID_t selfID );

  function AztecOO_StatusTestMaxIters_GetNumIters ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestMaxIters_GetNumIters')
    import :: c_int ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    integer(c_int)                                                            :: that
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name AztecOO_StatusTestResNorm interface
!! @{

  ! _________________ AztecOO_StatusTestResNorm interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTestResNorm_ID_t AztecOO_StatusTestResNorm_Degeneralize ( CTrilinos_Universal_ID_t id );

  function AztecOO_StatusTestResNorm_Degeneralize ( id ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_Degeneralize')
    import :: FT_AztecOO_StatusTestResNorm_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t)       ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t AztecOO_StatusTestResNorm_Generalize ( CT_AztecOO_StatusTestResNorm_ID_t id );

  function AztecOO_StatusTestResNorm_Generalize ( id ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                         :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusTestResNorm(const Epetra_Operator & Operator, const Epetra_Vector & LHS, 
  !!     const Epetra_Vector & RHS,double Tolerance);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusTestResNorm_ID_t AztecOO_StatusTestResNorm_Create ( CT_Epetra_Operator_ID_t OperatorID, 
  !!     CT_Epetra_Vector_ID_t LHSID, CT_Epetra_Vector_ID_t RHSID, double Tolerance );

  function AztecOO_StatusTestResNorm_Create ( OperatorID, LHSID, RHSID, Tolerance ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_Create')
    import :: FT_AztecOO_StatusTestResNorm_ID_t ,FT_Epetra_Operator_ID_t , &
          FT_Epetra_Vector_ID_t ,c_double
    
    type(FT_AztecOO_StatusTestResNorm_ID_t)                                  :: that
    type(FT_Epetra_Operator_ID_t)          ,intent(in)   ,value              :: OperatorID
    type(FT_Epetra_Vector_ID_t)            ,intent(in)   ,value              :: LHSID
    type(FT_Epetra_Vector_ID_t)            ,intent(in)   ,value              :: RHSID
    real(c_double)                         ,intent(in)   ,value              :: Tolerance
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~AztecOO_StatusTestResNorm();
  !> <BR> <BR> CTrilinos prototype:
  !! void AztecOO_StatusTestResNorm_Destroy ( CT_AztecOO_StatusTestResNorm_ID_t * selfID );

  subroutine AztecOO_StatusTestResNorm_Destroy ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_Destroy')
    import :: FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int DefineResForm( ResType TypeOfResidual, NormType TypeOfNorm, Epetra_Vector * Weights = 0);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_StatusTestResNorm_DefineResForm ( CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  !!     CT_ResType_E_t TypeOfResidual, CT_NormType_E_t TypeOfNorm, 
  !!     CT_Epetra_Vector_ID_t WeightsID );

  function AztecOO_StatusTestResNorm_DefineResForm ( selfID, TypeOfResidual, TypeOfNorm, &
        WeightsID ) result(that) bind(C,name='AztecOO_StatusTestResNorm_DefineResForm')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t ,FT_ResType_E_t ,FT_NormType_E_t , &
          FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                           :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    integer(FT_ResType_E_t)                ,intent(in)   ,value              :: TypeOfResidual
    integer(FT_NormType_E_t)               ,intent(in)   ,value              :: TypeOfNorm
    type(FT_Epetra_Vector_ID_t)            ,intent(in)   ,value              :: WeightsID
  end function


  !> <BR> Original C++ prototype:
  !! int DefineScaleForm( ScaleType TypeOfScaling, NormType TypeOfNorm, Epetra_Vector * Weights = 
  !!     0, double ScaleValue = 1.0);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_StatusTestResNorm_DefineScaleForm ( CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  !!     CT_ScaleType_E_t TypeOfScaling, CT_NormType_E_t TypeOfNorm, 
  !!     CT_Epetra_Vector_ID_t WeightsID, double ScaleValue );

  function AztecOO_StatusTestResNorm_DefineScaleForm ( selfID, TypeOfScaling, TypeOfNorm, &
        WeightsID, ScaleValue ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_DefineScaleForm')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t ,FT_ScaleType_E_t ,FT_NormType_E_t , &
          FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                           :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    integer(FT_ScaleType_E_t)              ,intent(in)   ,value              :: TypeOfScaling
    integer(FT_NormType_E_t)               ,intent(in)   ,value              :: TypeOfNorm
    type(FT_Epetra_Vector_ID_t)            ,intent(in)   ,value              :: WeightsID
    real(c_double)                         ,intent(in)   ,value              :: ScaleValue
  end function


  !> <BR> Original C++ prototype:
  !! int ResetTolerance(double Tolerance);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_StatusTestResNorm_ResetTolerance ( CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  !!     double Tolerance );

  function AztecOO_StatusTestResNorm_ResetTolerance ( selfID, Tolerance ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_ResetTolerance')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t ,c_double
    
    integer(c_int)                                                           :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                         ,intent(in)   ,value              :: Tolerance
  end function


  !> <BR> Original C++ prototype:
  !! int SetMaxNumExtraIterations(int maxNumExtraIterations);
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_StatusTestResNorm_SetMaxNumExtraIterations ( CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  !!     int maxNumExtraIterations );

  function AztecOO_StatusTestResNorm_SetMaxNumExtraIterations ( selfID, &
        maxNumExtraIterations ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_SetMaxNumExtraIterations')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t
    
    integer(c_int)                                                           :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                         ,intent(in)   ,value              :: maxNumExtraIterations
  end function


  !> <BR> Original C++ prototype:
  !! int GetMaxNumExtraIterations();
  !> <BR> <BR> CTrilinos prototype:
  !! int AztecOO_StatusTestResNorm_GetMaxNumExtraIterations ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  function AztecOO_StatusTestResNorm_GetMaxNumExtraIterations ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_GetMaxNumExtraIterations')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t
    
    integer(c_int)                                                           :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool ResidualVectorRequired() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean AztecOO_StatusTestResNorm_ResidualVectorRequired ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  function AztecOO_StatusTestResNorm_ResidualVectorRequired ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_ResidualVectorRequired')
    import :: FT_boolean_t ,FT_AztecOO_StatusTestResNorm_ID_t
    
    integer(FT_boolean_t)                                                    :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
  !!     double CurrentResNormEst, bool SolutionUpdated);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusType_E_t AztecOO_StatusTestResNorm_CheckStatus ( CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  !!     int CurrentIter, CT_Epetra_MultiVector_ID_t CurrentResVectorID, double CurrentResNormEst, 
  !!     boolean SolutionUpdated );

  function AztecOO_StatusTestResNorm_CheckStatus ( selfID, CurrentIter, CurrentResVectorID, &
        CurrentResNormEst, SolutionUpdated ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_CheckStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestResNorm_ID_t ,c_int , &
          FT_Epetra_MultiVector_ID_t ,c_double ,FT_boolean_t
    
    integer(FT_AztecOO_StatusType_E_t)                                       :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                         ,intent(in)   ,value              :: CurrentIter
    type(FT_Epetra_MultiVector_ID_t)       ,intent(in)   ,value              :: CurrentResVectorID
    real(c_double)                         ,intent(in)   ,value              :: CurrentResNormEst
    integer(FT_boolean_t)                  ,intent(in)   ,value              :: SolutionUpdated
  end function


  !> <BR> Original C++ prototype:
  !! AztecOO_StatusType GetStatus() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_AztecOO_StatusType_E_t AztecOO_StatusTestResNorm_GetStatus ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  function AztecOO_StatusTestResNorm_GetStatus ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_GetStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestResNorm_ID_t
    
    integer(FT_AztecOO_StatusType_E_t)                                       :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void ResetStatus();
  !> <BR> <BR> CTrilinos prototype:
  !! void AztecOO_StatusTestResNorm_ResetStatus ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  subroutine AztecOO_StatusTestResNorm_ResetStatus ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_ResetStatus')
    import :: FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! double GetTolerance() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double AztecOO_StatusTestResNorm_GetTolerance ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  function AztecOO_StatusTestResNorm_GetTolerance ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_GetTolerance')
    import :: c_double ,FT_AztecOO_StatusTestResNorm_ID_t
    
    real(c_double)                                                           :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double GetTestValue() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double AztecOO_StatusTestResNorm_GetTestValue ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  function AztecOO_StatusTestResNorm_GetTestValue ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_GetTestValue')
    import :: c_double ,FT_AztecOO_StatusTestResNorm_ID_t
    
    real(c_double)                                                           :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double GetResNormValue() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double AztecOO_StatusTestResNorm_GetResNormValue ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  function AztecOO_StatusTestResNorm_GetResNormValue ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_GetResNormValue')
    import :: c_double ,FT_AztecOO_StatusTestResNorm_ID_t
    
    real(c_double)                                                           :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double GetScaledNormValue() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double AztecOO_StatusTestResNorm_GetScaledNormValue ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  function AztecOO_StatusTestResNorm_GetScaledNormValue ( selfID ) result(that) &
        bind(C,name='AztecOO_StatusTestResNorm_GetScaledNormValue')
    import :: c_double ,FT_AztecOO_StatusTestResNorm_ID_t
    
    real(c_double)                                                           :: that
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


  end interface
end module foraztecoo

#endif /* HAVE_FORTRILINOS_AZTECOO */
