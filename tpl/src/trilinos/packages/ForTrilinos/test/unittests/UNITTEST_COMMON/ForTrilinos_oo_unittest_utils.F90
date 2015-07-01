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

module ForTrilinos_oo_unittest_utils
#include "ForTrilinos_config.h"
  use iso_c_binding ,only : c_int        ! Kind parameter (precision specifier)
  use ForTrilinos_enums
#ifdef HAVE_MPI
  use mpi, only:MPI_COMM_WORLD
  use FEpetra_MpiComm,only:Epetra_MpiComm
#else 
  use FEpetra_SerialComm,only:Epetra_SerialComm
#endif
  implicit none                          ! Prevent implicit typing
!  public :: comm

!#ifdef HAVE_MPI
!  type(Epetra_MpiComm) :: comm
!#else
!   type(Epetra_SerialComm) :: comm
!#endif
!  end type

  contains

#ifdef HAVE_MPI

  ! /*! Create an Epetra_MpiComm */
  type(Epetra_MpiComm) function UnitTest_EpetraComm_Create() 
!    use mpi,only:MPI_COMM_WORLD
    UnitTest_EpetraComm_Create = Epetra_MpiComm(MPI_COMM_WORLD)
  end function
#else
  ! /*! Create an Epetra_SerialComm */
  type(Epetra_SerialComm) function UnitTest_EpetraComm_Create() 
    UnitTest_EpetraComm_Create = Epetra_SerialComm()
  end function

#endif

end module ForTrilinos_oo_unittest_utils
