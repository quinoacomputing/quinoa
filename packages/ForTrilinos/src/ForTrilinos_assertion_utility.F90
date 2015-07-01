!*********************************************************************
! ForTrilinos: Object-Oriented Fortran 2003 interface to Trilinos
!		 Copyright 2010 Sandia Corporation 
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
! Questions? Contact Karla Morris  (knmorri@sandia.gov)
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************

module ForTrilinos_assertion_utility
#include "ForTrilinos_config.h"
  use iso_fortran_env ,only : error_unit  
  implicit none
  private
  public :: error_message,assert,assert_identical

#ifdef ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS
  integer ,parameter :: max_string_length=256
#endif /* ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS */

  type error_message
    private
#ifdef ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS
    character(len=max_string_length) :: string ! gfortran 4.7.0 workaround
#else
    character(:) ,allocatable :: string
#endif /* ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS */
  contains 
    procedure :: text
  end type

  interface error_message ! constructor
    module procedure new_message
  end interface 

  interface assert
    module procedure scalar_assert,vector_assert
  end interface

contains

  type(error_message) function new_message(message)
    character(len=*), intent(in) :: message
    new_message%string = message 
  end function

  function text(this)
    class(error_message) ,intent(in) :: this
#ifdef ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS
    character(len=max_string_length) :: text
#else
    character(:) ,allocatable :: text
    if (allocated(this%string)) then
#endif /* ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS */
       text = this%string
#ifndef ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS
    else
       text = 'No error message provided.'
    end if
#endif /* ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS */
  end function

  subroutine scalar_assert(assertion,text)
    logical ,intent(in) :: assertion
    type(error_message) ,intent(in) :: text
    if (.not. assertion) then
      write(error_unit,fmt='(31a)',advance="no") 'Assertion failed with message: '
#ifndef ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS
      if (allocated(text%string)) then
#endif
        write(error_unit,*) text%string
#ifndef ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS
      else
        write(error_unit,*) '(no message provided).'
      end if
#endif
    end if
  end subroutine

  subroutine vector_assert(assertion,text)
    logical ,dimension(:) ,intent(in) :: assertion
    type(error_message) ,dimension(:) ,intent(in) :: text
    integer :: i
    logical :: any_failures 
    call assert_identical( [size(assertion),size(text)] )
    any_failures=.false.
    do i=1,size(assertion)
      if (.not. assertion(i)) then
        any_failures=.true.
        write(error_unit,fmt='(31a)',advance="no") 'Assertion failed with message: '
       !if (allocated(text(i)%string)) then
          write(error_unit,*) text(i)%string
       !else
       !  write(error_unit,*) '(no message provided).'
       !end if
      end if
    end do
    if (any_failures) stop 'Execution halted on failed assertion(s)!'
  end subroutine

  subroutine assert_identical(integers)
    integer ,dimension(:) ,intent(in) :: integers
    integer :: i
    logical :: any_mismatches
    any_mismatches = .false.
    do i=2,size(integers)
      if (integers(i) /= integers(1)) then
        any_mismatches = .true.
        write(error_unit,*) &
        'Value ',i,' does not match expected value ',integers(1)
      end if
    end do
    if (any_mismatches) stop 'Execution halted on failed assertion!'
  end subroutine
end module
