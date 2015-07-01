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

module ForTrilinos_table_man
#include "ForTrilinos_config.h"
  use iso_c_binding ,only : c_int        ! Kind parameter (precision specifier)
  use ForTrilinos_enums

  implicit none                          ! Prevent implicit typing

  interface

  ! /*! Copies the RCP from one table into a second table. The new ID
  !  *  will be returned from the function. Both the old and the new
  !  *  IDs will need to be removed from the tables in order to destroy
  !  *  the object. */
  ! CTrilinos_Universal_ID_t CT_Alias(CTrilinos_Universal_ID_t aid,
  ! CTrilinos_Table_ID_t new_table);

  type(ForTrilinos_Universal_ID_t) function CT_Alias( selfID, new_table ) &
        bind(C,name='CT_Alias')
    import :: ForTrilinos_Universal_ID_t, ForTrilinos_Table_ID_t

    type(ForTrilinos_Universal_ID_t),intent(in),value :: selfID
    integer(ForTrilinos_Table_ID_t) ,intent(in),value :: new_table
  end function

  ! /*! Removes the RCP from one table and puts it in another. *aid will
  !  *  hold the new struct value afterward. Only the new RCP will need
  !  *  to be removed in order to destroy the object. */
  ! void CT_Migrate(CTrilinos_Universal_ID_t *aid, CTrilinos_Table_ID_t new_table);

  subroutine CT_Migrate( selfID, new_table ) &
        bind(C,name='CT_Migrate')
    import :: ForTrilinos_Universal_ID_t, ForTrilinos_Table_ID_t

    type(ForTrilinos_Universal_ID_t)                  :: selfID
    integer(ForTrilinos_Table_ID_t) ,intent(in),value :: new_table
  end subroutine

  ! /*! Checks to see if the underlying object referenced by a table
  !  *  entry is dynamic_cast'able to a given type (can be used to
  !  *  distinguish, e.g., an Epetra_SerialComm from an Epetra_MpiComm). */
  ! boolean CT_TypeCheck(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t type);

  subroutine CT_TypeCheck( selfID, typeid ) &
        bind(C,name='CT_TypeCheck')
    import :: ForTrilinos_Universal_ID_t, ForTrilinos_Table_ID_t

    type(ForTrilinos_Universal_ID_t),intent(in),value :: selfID
    integer(ForTrilinos_Table_ID_t) ,intent(in),value :: typeid
  end subroutine

  end interface
end module ForTrilinos_table_man
