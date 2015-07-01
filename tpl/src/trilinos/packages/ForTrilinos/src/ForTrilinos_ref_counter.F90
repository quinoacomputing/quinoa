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
! Questions? Contact Karla Morris  (knmorri@sandia.gov)
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************

module ForTrilinos_ref_counter
#include "ForTrilinos_config.h"
  use ForTrilinos_hermetic ,only : hermetic
  use ForTrilinos_assertion_utility ,only : assert,error_message
  use ForTrilinos_error ,only : error
  implicit none
  private
  public :: ref_counter
  type ref_counter
      private
      integer, pointer :: count => null()
      class(hermetic), pointer :: obj => null()
  contains
      procedure, non_overridable :: grab
      procedure, non_overridable :: release
      procedure :: assign
#ifndef ForTrilinos_DISABLE_FINAL_SUBROUTINES
      final     :: finalize_ref_counter
#endif /* ForTrilinos_DISABLE_FINAL_SUBROUTINES */
      generic   :: assignment(=) => assign
  end type

  interface ref_counter
      module procedure new_ref_counter
  end interface

contains

  subroutine grab(this)
    class(ref_counter), intent(inout) :: this
    call assert( associated(this%count), error_message('Ref_counter%grab: count not associated.') )
!!$    call assert( [associated(this%count)], [error_message('Ref_counter%grab: count not associated.')] )
    this%count = this%count + 1
  end subroutine

  subroutine release(this)
    class (ref_counter), intent(inout) :: this
    integer  :: status
    type(error) :: ierr
    call assert( associated(this%count), error_message('Ref_counter%release: count not associated.') )
    call assert( (this%count>=1), error_message('Ref_counter%release: non-positive count.') )
!!$    call assert( [associated(this%count)], [error_message('Ref_counter%release: count not associated.')] )
!!$    call assert( [this%count>=1], [error_message('Ref_counter%release: non-positive count.')] )
    this%count = this%count - 1
    if (this%count == 0) then
      call this%obj%ctrilinos_delete
      deallocate(this%count,stat=status)
      ierr=error(status,'Ref_counter%release: this%count deallocation failed.')
      call ierr%check_success()
      !This would be use once allocate(this%obj,mold=object) is working
      !deallocate(this%obj,stat=status)
      !ierr=error(status,'Ref_counter%release: this%obj deallocation failed.')
      !call ierr%check_success()
      !This would be use once allocate(this%obj,mold=object) is working
      this%obj => null() ! workaround for lack of mold feature deals with invalid value in this%count of a reference counted component.
    else
      this%count=>null()
      this%obj=>null()
    end if
  end subroutine

  subroutine assign (lhs, rhs)
    class (ref_counter), intent(inout) :: lhs
    class (ref_counter), intent(in) :: rhs
    call assert( associated(rhs%count), error_message('Ref_counter%assign: rhs%count not associated.') )
    call assert( associated(rhs%obj), error_message('Ref_counter%assign: rhs%obj not associated.') )
    lhs%count => rhs%count
    lhs%obj => rhs%obj
    if (associated(lhs%count)) call lhs%grab
  end subroutine

  recursive subroutine finalize_ref_counter (this)
    type(ref_counter), intent(inout) :: this
    if (associated(this%count)) call this%release
  end subroutine

  function new_ref_counter(object)
    class(hermetic), intent(in) :: object
    integer :: status
    type(error) :: ierr
    type(ref_counter), allocatable :: new_ref_counter

    allocate(new_ref_counter,stat=status)
    ierr=error(status,'new_ref_counter: new_ref_counter allocation failed.')
    call ierr%check_success()

    allocate(new_ref_counter%count, source=0,stat=status)
    ierr=error(status,'new_ref_counter: count allocation failed.')
    call ierr%check_success()

    allocate(new_ref_counter%obj,source=object,stat=status)
    ierr=error(status,'new_ref_counter: obj allocation failed.')
    call ierr%check_success()

    call new_ref_counter%grab
  end function
end module
