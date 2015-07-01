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
module ForTrilinos_enum_wrappers
  use iso_c_binding ,only : c_int        ! Kind parameter (precision specifier)
  implicit none                          ! Prevent implicit typing

  ! Epetra_DataAccess
  integer(kind(c_int)) ,parameter :: FT_Epetra_DataAccess_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_Epetra_DataAccess_E_Copy,                                    &
      FT_Epetra_DataAccess_E_View                                     
  end enum

  ! Epetra_CombineMode
  integer(kind(c_int)) ,parameter :: FT_Epetra_CombineMode_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_Epetra_CombineMode_E_Add,                                    &
      FT_Epetra_CombineMode_E_Zero,                                   &
      FT_Epetra_CombineMode_E_Insert,                                 &
      FT_Epetra_CombineMode_E_InsertAdd,                              &
      FT_Epetra_CombineMode_E_Average,                                &
      FT_Epetra_CombineMode_E_AbsMax                                  
  end enum

  ! ProblemDifficultyLevel
  integer(kind(c_int)) ,parameter :: FT_ProblemDifficultyLevel_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_ProblemDifficultyLevel_E_easy,                               &
      FT_ProblemDifficultyLevel_E_moderate,                           &
      FT_ProblemDifficultyLevel_E_hard,                               &
      FT_ProblemDifficultyLevel_E_unsure                              
  end enum

  ! Teuchos::EValidateUsed
  integer(kind(c_int)) ,parameter :: FT_EValidateUsed_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_EValidateUsed_E_VALIDATE_USED_ENABLED,                       &
      FT_EValidateUsed_E_VALIDATE_USED_DISABLED                       
  end enum

  ! Teuchos::EValidateDefaults
  integer(kind(c_int)) ,parameter :: FT_EValidateDefaults_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_EValidateDefaults_E_VALIDATE_DEFAULTS_ENABLED,               &
      FT_EValidateDefaults_E_VALIDATE_DEFAULTS_DISABLED               
  end enum

  ! Teuchos::CommandLineProcessor::EParseCommandLineReturn
  integer(kind(c_int)) ,parameter :: FT_EParseCommandLineReturn_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_EParseCommandLineReturn_E_PARSE_SUCCESSFUL = 0,              &
      FT_EParseCommandLineReturn_E_PARSE_HELP_PRINTED = 1,            &
      FT_EParseCommandLineReturn_E_PARSE_UNRECOGNIZED_OPTION = 2      
  end enum

  ! Teuchos::CommandLineProcessor::EOptType
  integer(kind(c_int)) ,parameter :: FT_EOptType_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_EOptType_E_OPT_NONE,                                         &
      FT_EOptType_E_OPT_BOOL_TRUE,                                    &
      FT_EOptType_E_OPT_BOOL_FALSE,                                   &
      FT_EOptType_E_OPT_INT,                                          &
      FT_EOptType_E_OPT_DOUBLE,                                       &
      FT_EOptType_E_OPT_STRING,                                       &
      FT_EOptType_E_OPT_ENUM_INT                                      
  end enum

#ifdef HAVE_FORTRILINOS_AZTECOO
  ! AztecOO_StatusType
  integer(kind(c_int)) ,parameter :: FT_AztecOO_StatusType_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_AztecOO_StatusType_E_Unchecked = 2,                          &
      FT_AztecOO_StatusType_E_Unconverged = 1,                        &
      FT_AztecOO_StatusType_E_Converged = 0,                          &
      FT_AztecOO_StatusType_E_Failed = -1,                            &
      FT_AztecOO_StatusType_E_NaN = -2,                               &
      FT_AztecOO_StatusType_E_PartialFailed = -3                      
  end enum
#endif /* HAVE_FORTRILINOS_AZTECOO */

#ifdef HAVE_FORTRILINOS_AZTECOO
  ! AztecOO_StatusTestCombo::ComboType
  integer(kind(c_int)) ,parameter :: FT_ComboType_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_ComboType_E_AND,                                             &
      FT_ComboType_E_OR,                                              &
      FT_ComboType_E_SEQ                                              
  end enum
#endif /* HAVE_FORTRILINOS_AZTECOO */

#ifdef HAVE_FORTRILINOS_AZTECOO
  ! AztecOO_StatusTestResNorm::ResType
  integer(kind(c_int)) ,parameter :: FT_ResType_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_ResType_E_Implicit,                                          &
      FT_ResType_E_Explicit                                           
  end enum
#endif /* HAVE_FORTRILINOS_AZTECOO */

#ifdef HAVE_FORTRILINOS_AZTECOO
  ! AztecOO_StatusTestResNorm::NormType
  integer(kind(c_int)) ,parameter :: FT_NormType_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_NormType_E_OneNorm,                                          &
      FT_NormType_E_TwoNorm,                                          &
      FT_NormType_E_InfNorm                                           
  end enum
#endif /* HAVE_FORTRILINOS_AZTECOO */

#ifdef HAVE_FORTRILINOS_AZTECOO
  ! AztecOO_StatusTestResNorm::ScaleType
  integer(kind(c_int)) ,parameter :: FT_ScaleType_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_ScaleType_E_NormOfRHS,                                       &
      FT_ScaleType_E_NormOfInitRes,                                   &
      FT_ScaleType_E_None,                                            &
      FT_ScaleType_E_UserProvided                                     
  end enum
#endif /* HAVE_FORTRILINOS_AZTECOO */

#ifdef HAVE_FORTRILINOS_IFPACK
  ! Ifpack::EPrecType
  integer(kind(c_int)) ,parameter :: FT_EPrecType_E_t = c_int

  enum ,bind(C)
    enumerator ::                                                     &
      FT_EPrecType_E_POINT_RELAXATION,                                &
      FT_EPrecType_E_POINT_RELAXATION_STAND_ALONE,                    &
      FT_EPrecType_E_BLOCK_RELAXATION,                                &
      FT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE,                    &
      FT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_ILU,                &
      FT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_AMESOS,             &
      FT_EPrecType_E_BLOCK_RELAXATION_AMESOS,                         &
      FT_EPrecType_E_AMESOS,                                          &
      FT_EPrecType_E_AMESOS_STAND_ALONE,                              &
      FT_EPrecType_E_IC,                                              &
      FT_EPrecType_E_IC_STAND_ALONE,                                  &
      FT_EPrecType_E_ICT,                                             &
      FT_EPrecType_E_ICT_STAND_ALONE,                                 &
      FT_EPrecType_E_ILU,                                             &
      FT_EPrecType_E_ILU_STAND_ALONE,                                 &
      FT_EPrecType_E_ILUT,                                            &
      FT_EPrecType_E_ILUT_STAND_ALONE,                                &
      FT_EPrecType_E_SPARSKIT,                                        &
      FT_EPrecType_E_HIPS,                                            &
      FT_EPrecType_E_HYPRE,                                           &
      FT_EPrecType_E_CHEBYSHEV                                        
  end enum
#endif /* HAVE_FORTRILINOS_IFPACK */

end module ForTrilinos_enum_wrappers
