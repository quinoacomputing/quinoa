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
#ifdef HAVE_FORTRILINOS_AMESOS

module foramesos
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  use ForTrilinos_enum_wrappers
  implicit none   ! Prevent implicit typing
#ifdef HAVE_MPI
#include <mpif.h>
#endif

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/amesos/CAmesos*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

!> @name Amesos_BaseSolver interface
!! @{

  ! _________________ Amesos_BaseSolver interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Amesos_BaseSolver_ID_t Amesos_BaseSolver_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Amesos_BaseSolver_Degeneralize ( id ) result(that) &
        bind(C,name='Amesos_BaseSolver_Degeneralize')
    import :: FT_Amesos_BaseSolver_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Amesos_BaseSolver_Generalize ( CT_Amesos_BaseSolver_ID_t id );

  function Amesos_BaseSolver_Generalize ( id ) result(that) &
        bind(C,name='Amesos_BaseSolver_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Amesos_BaseSolver();
  !> <BR> <BR> CTrilinos prototype:
  !! void Amesos_BaseSolver_Destroy ( CT_Amesos_BaseSolver_ID_t * selfID );

  subroutine Amesos_BaseSolver_Destroy ( selfID ) bind(C,name='Amesos_BaseSolver_Destroy')
    import :: FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual int SymbolicFactorization() = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Amesos_BaseSolver_SymbolicFactorization ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_SymbolicFactorization ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_SymbolicFactorization')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumericFactorization() = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Amesos_BaseSolver_NumericFactorization ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_NumericFactorization ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_NumericFactorization')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Solve() = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Amesos_BaseSolver_Solve ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_Solve ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_Solve')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int SetUseTranspose(bool UseTranspose) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Amesos_BaseSolver_SetUseTranspose ( CT_Amesos_BaseSolver_ID_t selfID, 
  !!     boolean UseTranspose );

  function Amesos_BaseSolver_SetUseTranspose ( selfID, UseTranspose ) result(that) &
        bind(C,name='Amesos_BaseSolver_SetUseTranspose')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t ,FT_boolean_t
    
    integer(c_int)                                                   :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)          ,intent(in)   ,value              :: UseTranspose
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool UseTranspose() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Amesos_BaseSolver_UseTranspose ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_UseTranspose ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_UseTranspose')
    import :: FT_boolean_t ,FT_Amesos_BaseSolver_ID_t
    
    integer(FT_boolean_t)                                            :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int SetParameters( Teuchos::ParameterList &ParameterList ) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Amesos_BaseSolver_SetParameters ( CT_Amesos_BaseSolver_ID_t selfID, 
  !!     CT_Teuchos_ParameterList_ID_t ParameterListID );

  function Amesos_BaseSolver_SetParameters ( selfID, ParameterListID ) result(that) &
        bind(C,name='Amesos_BaseSolver_SetParameters')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t ,FT_Teuchos_ParameterList_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ParameterListID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_LinearProblem* GetProblem() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_LinearProblem_ID_t Amesos_BaseSolver_GetProblem ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_GetProblem ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_GetProblem')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t)                                  :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool MatrixShapeOK() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Amesos_BaseSolver_MatrixShapeOK ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_MatrixShapeOK ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_MatrixShapeOK')
    import :: FT_boolean_t ,FT_Amesos_BaseSolver_ID_t
    
    integer(FT_boolean_t)                                            :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_Comm & Comm() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Amesos_BaseSolver_Comm ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_Comm ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_Comm')
    import :: FT_Epetra_Comm_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                        :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumSymbolicFact() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Amesos_BaseSolver_NumSymbolicFact ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_NumSymbolicFact ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_NumSymbolicFact')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumNumericFact() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Amesos_BaseSolver_NumNumericFact ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_NumNumericFact ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_NumNumericFact')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumSolve() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Amesos_BaseSolver_NumSolve ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_NumSolve ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_NumSolve')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual void PrintStatus() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! void Amesos_BaseSolver_PrintStatus ( CT_Amesos_BaseSolver_ID_t selfID );

  subroutine Amesos_BaseSolver_PrintStatus ( selfID ) &
        bind(C,name='Amesos_BaseSolver_PrintStatus')
    import :: FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual void PrintTiming() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! void Amesos_BaseSolver_PrintTiming ( CT_Amesos_BaseSolver_ID_t selfID );

  subroutine Amesos_BaseSolver_PrintTiming ( selfID ) &
        bind(C,name='Amesos_BaseSolver_PrintTiming')
    import :: FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  !> <BR> <BR> CTrilinos prototype:
  !! void Amesos_BaseSolver_setParameterList ( CT_Amesos_BaseSolver_ID_t selfID, 
  !!     CT_Teuchos_ParameterList_ID_t paramListID );

  subroutine Amesos_BaseSolver_setParameterList ( selfID, paramListID ) &
        bind(C,name='Amesos_BaseSolver_setParameterList')
    import :: FT_Amesos_BaseSolver_ID_t ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: paramListID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_getNonconstParameterList ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_getNonconstParameterList ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_getNonconstParameterList')
    import :: FT_Teuchos_ParameterList_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_unsetParameterList ( CT_Amesos_BaseSolver_ID_t selfID );

  function Amesos_BaseSolver_unsetParameterList ( selfID ) result(that) &
        bind(C,name='Amesos_BaseSolver_unsetParameterList')
    import :: FT_Teuchos_ParameterList_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual void GetTiming( Teuchos::ParameterList &TimingParameterList ) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Amesos_BaseSolver_GetTiming ( CT_Amesos_BaseSolver_ID_t selfID, 
  !!     CT_Teuchos_ParameterList_ID_t TimingParameterListID );

  subroutine Amesos_BaseSolver_GetTiming ( selfID, TimingParameterListID ) &
        bind(C,name='Amesos_BaseSolver_GetTiming')
    import :: FT_Amesos_BaseSolver_ID_t ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: TimingParameterListID
  end subroutine


!> @}


!> @name Amesos interface
!! @{

  ! _________________ Amesos interface bodies _________________


  !> <BR> Original C++ prototype:
  !! Amesos();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Amesos_ID_t Amesos_Create (  );

  function Amesos_Create (  ) result(that) bind(C,name='Amesos_Create')
    import :: FT_Amesos_ID_t
    
    type(FT_Amesos_ID_t)                                          :: that
  end function


  !> <BR> Original C++ prototype:
  !! ~Amesos();
  !> <BR> <BR> CTrilinos prototype:
  !! void Amesos_Destroy ( CT_Amesos_ID_t * selfID );

  subroutine Amesos_Destroy ( selfID ) bind(C,name='Amesos_Destroy')
    import :: FT_Amesos_ID_t
    
    type(FT_Amesos_ID_t)                                          :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! Amesos_BaseSolver *Create(const char *ClassType, const Epetra_LinearProblem& LinearProblem );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Amesos_BaseSolver_ID_t Amesos_CreateSolver ( CT_Amesos_ID_t selfID, const char * ClassType, 
  !!     CT_Epetra_LinearProblem_ID_t LinearProblemID );

  function Amesos_CreateSolver ( selfID, ClassType, LinearProblemID ) result(that) &
        bind(C,name='Amesos_CreateSolver')
    import :: FT_Amesos_BaseSolver_ID_t ,FT_Amesos_ID_t ,c_char , &
          FT_Epetra_LinearProblem_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t)                                  :: that
    type(FT_Amesos_ID_t)        ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: ClassType
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: LinearProblemID
  end function


  !> <BR> Original C++ prototype:
  !! bool Query(const char * ClassType);
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Amesos_Query ( CT_Amesos_ID_t selfID, const char * ClassType );

  function Amesos_Query ( selfID, ClassType ) result(that) bind(C,name='Amesos_Query')
    import :: FT_boolean_t ,FT_Amesos_ID_t ,c_char
    
    integer(FT_boolean_t)                                         :: that
    type(FT_Amesos_ID_t)        ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: ClassType
  end function


  !> <BR> Original C++ prototype:
  !! static Teuchos::ParameterList GetValidParameters();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Amesos_GetValidParameters (  );

  function Amesos_GetValidParameters (  ) result(that) &
        bind(C,name='Amesos_GetValidParameters')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
  end function


!> @}


  end interface
end module foramesos

#endif /* HAVE_FORTRILINOS_AMESOS */
