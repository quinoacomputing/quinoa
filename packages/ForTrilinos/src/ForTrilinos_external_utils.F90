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

module ForTrilinos_external_utils
#include "ForTrilinos_config.h"
  use iso_c_binding ,only : c_int        ! Kind parameter (precision specifier)
  use ForTrilinos_enums

  implicit none                          ! Prevent implicit typing

  interface

#ifdef HAVE_MPI

  ! /*! Create an Epetra_MpiComm from Fortran */
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Fortran_Create ( int fcomm );

  type(FT_Epetra_MpiComm_ID_t) function Epetra_MpiComm_Fortran_Create( fcomm ) &
        bind(C,name='Epetra_MpiComm_Fortran_Create')
    import :: FT_Epetra_MpiComm_ID_t, c_int

    integer(c_int),intent(in),value :: fcomm
  end function

#endif

  subroutine ForTrilinos_CleanSlate( ) &
        bind(C,name='CTrilinos_CleanSlate')
  end subroutine

  end interface
end module ForTrilinos_external_utils
