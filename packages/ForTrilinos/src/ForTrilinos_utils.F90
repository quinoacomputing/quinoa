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


module fortrilinos_utils
#include "ForTrilinos_config.h"
  implicit none
  private
  public :: c_strlen
  public :: fortran_string
  public :: generalize_all
  public :: valid_kind_parameters
contains
  
  ! Fortran-style strings take two forms: one with no C counterpart and a second that is
  ! an array of character values.  A C string is interoperable with the second type of 
  ! Fortran string with the primary distinction being that Fortran strings carry their
  ! own string-length information, whereas the length of C strings is indicated by a 
  ! terminating null character. This function takes in a c_string and returns the 
  ! equivalent Fortran string.

  function fortran_string(c_string)
    use ,intrinsic :: iso_c_binding ,only: c_char,c_null_char
    character(kind=c_char) ,intent(in) :: c_string(*)
    character(:) ,allocatable :: fortran_string
    integer :: i
    i=1
    fortran_string = c_string(i)
    do while(c_string(i)/=c_null_char)
      i=i+1
      fortran_string = fortran_string // c_string(i)
    end do
  end function

  ! This is a Fortran implementation of the C strlen function that returns
  ! the C string length by searching for the null terminator.  The length
  ! differs from the result of the Fortran len intrinsic by a unit decrement
  ! corresponding to the position of the terminating null character.

  function c_strlen(c_string)
    use ,intrinsic :: iso_c_binding ,only: c_char ,c_int  ,c_null_char
    character(kind=c_char) ,intent(in) :: c_string(*)
    integer(c_int) :: c_strlen
    c_strlen = len(fortran_string(c_string)) - len(c_null_char)
  end function

  ! This procedure converts a C string to a variable-length Fortran string.


  ! This is a Fortran implementation of the functionality in the Epetra_*_Abstract procedures
  ! in each CTrilinos/src/CEpetra* file.  It effectively casts any given Epetra derived type
  ! to a general type that can represent any of the Epetra derived types in 
  ! ForTrilinos/src/ForTrilinos_enums.F90.

  type(ForTrilinos_Universal_ID_t) function generalize_all(object_id) bind(C) 
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr) ,value :: object_id
    type(ForTrilinos_Universal_ID_t), pointer :: local_ptr
    call c_f_pointer (object_id, local_ptr)
    generalize_all = local_ptr
  end function

  ! This procedure checks the values of parameters required to interoperate with CTrilinos.  The Fortran 2003 
  ! standard requires that these parameters be defined in the intrinsic module iso_c_binding with values
  ! that communicate meanings specified in the standard.  This procedure's quoted interpretations of these 
  ! values are largely excerpted from the standard.  This procedure returns true if all of the interoperating
  ! Fortran kind  parameters required by ForTrilinos have a corresponding C type defined by the companion 
  ! C processor.  Otherwise, return false.  (For purposes of ForTrilinos, the Fortran standard's use of the
  ! word 'processor' is interpreted as denoting the combination of a compiler, an operating system, and a
  ! hardware architecture.)

  logical function valid_kind_parameters(verbose)
    use ,intrinsic :: iso_c_binding ,only : &
      c_int,c_char,c_double,c_ptr,c_long,c_bool,c_null_char
    use ,intrinsic :: iso_fortran_env ,only : error_unit ,output_unit
    logical ,optional :: verbose
    logical           :: verbose_output


    character(len=*) ,parameter :: no_fortran_kind= &
      'The companion C processor defines the corresponding C type, &
      & but there is no interoperating Fortran processor kind.'
    character(len=*) ,parameter :: no_c_type= &
      'The C processor does not define the corresponding C type.'
    character(len=*) ,parameter :: interoperable = &
      'This interoperating Fortran kind has a corresponding C type.'
    character(len=*) ,parameter :: imprecise= &
      'The C processor type does not have a precision equal to the precision of any of the Fortran processor real kinds.'
    character(len=*) ,parameter :: limited= &
      'The C processor type does not have a range equal to the range of any of the Fortran processor real kinds.'
    character(len=*) ,parameter :: limited_and_imprecise= &
      'The C processor type has neither the precision nor range of any of the Fortran processor real kinds.'
    character(len=*) ,parameter :: not_interoperable_nonspecific = &
      'There is no interoperating Fortran processor kind for unspecified reasons.'

    valid_kind_parameters = .true. ! default return value

    if (present(verbose)) then 
      verbose_output=verbose
    else
      verbose_output=.false.
    end if

    select case(c_long)
      case(-1)
        write(error_unit ,fmt='(2a)') 'c_long error: ',no_fortran_kind
        valid_kind_parameters = .false.
      case(-2)
        write(error_unit ,fmt='(2a)') 'c_long error: ',no_c_type
        valid_kind_parameters = .false.
      case default
        if (verbose_output) write(output_unit,fmt='(2a)') 'c_long: ',interoperable
    end select

    select case(c_int)
      case(-1)
        write(error_unit ,fmt='(2a)') 'c_int error: ',no_fortran_kind
        valid_kind_parameters = .false.
      case(-2)
        write(error_unit ,fmt='(2a)') 'c_int error: ',no_c_type
        valid_kind_parameters = .false.
      case default
        if (verbose_output) write(output_unit,fmt='(2a)') 'c_int: ',interoperable
    end select
 
    select case(c_double)
      case(-1)
        write(error_unit ,fmt='(2a)') 'c_double error: ',imprecise
        valid_kind_parameters = .false.
      case(-2)
        write(error_unit ,fmt='(2a)') 'c_double error: ',limited
        valid_kind_parameters = .false.
      case(-3)
        write(error_unit ,fmt='(2a)') 'c_double error: ',limited_and_imprecise
        valid_kind_parameters = .false.
      case(-4)
        write(error_unit ,fmt='(2a)') 'c_double error: ',not_interoperable_nonspecific
        valid_kind_parameters = .false.
      case default
        if (verbose_output) write(output_unit,fmt='(2a)') 'c_double: ',interoperable
    end select
 
    select case(c_bool)
      case(-1)
        write(error_unit ,fmt='(a)') 'c_bool error: invalid value for a logical kind parameter on the processor.'
        valid_kind_parameters = .false.
      case default
        if (verbose_output) write(output_unit ,fmt='(a)') 'c_bool:  valid value for a logical kind parameter on the processor.'
    end select
  
    select case(c_char)
      case(-1)
        write(error_unit ,fmt='(a)') 'c_char error: invalid value for a character kind parameter on the processor.'
        valid_kind_parameters = .false.
      case default
        if (verbose_output) write(output_unit ,fmt='(a)') 'c_char:  valid value for a character kind parameter on the processor.'
    end select
  end function

end module fortrilinos_utils
