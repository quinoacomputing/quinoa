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


module FEpetra_SrcDistObject
  use ForTrilinos_enums !,only: FT_Epetra_SrcDistObject_ID_t,ForTrilinos_Universal_ID_t,FT_Epetra_BlockMap_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only : universal
  use ForTrilinos_error
  use FEpetra_BlockMap, only: Epetra_BlockMap
  use forepetra
  implicit none
  !private               ! Hide everything by default
  !public :: Epetra_SrcDistObject ! Expose type/methods
 
  !> <BR> Epetra_SrcDistObject: A class for supporting flexible source distributed objects for import/export operations.
  !> @{

  !> @brief The Epetra_SrcDistObject is a base class for all Epetra distributed global objects that are potential source objects for the general Epetra_DistObject class.  
  !!It provides a way to send a very general distributed object as the potential source object for an import or export object.  

  type ,abstract :: Epetra_SrcDistObject !,extends(universal) :: Epetra_SrcDistObject
  contains
    ! Developers only
    !procedure(map_interface),deferred  :: map
  end type
  
  
  !abstract interface
  ! We need a CTrilinos function that tells us the type of the underlying object
  !function map_interface(this) result(map)  
  !  import :: Epetra_BlockMap, Epetra_SrcDistObject
  !  class(Epetra_BlockMap), allocatable :: map  
  !  class(Epetra_SrcDistObject), intent(in) :: this
  !end function
  !end interface
 
end module 
