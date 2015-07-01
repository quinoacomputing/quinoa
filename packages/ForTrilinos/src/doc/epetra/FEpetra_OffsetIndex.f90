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


module FEpetra_OffsetIndex
  use ForTrilinos_enums ,only: FT_Epetra_OffsetIndex_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_Import ,only: Epetra_Import
  use FEpetra_Export ,only: Epetra_Export
  use FEpetra_CrsGraph ,only: Epetra_CrsGraph
  use iso_c_binding     ,only: c_int,c_double,c_char
  use forepetra
  implicit none
  !private                      ! Hide everything by default
  !public :: Epetra_OffsetIndex ! Expose type/constructors/methods

  !> <BR> Epetra_OffsetIndex: This class builds index for efficient mapping of data from one Epetra_CrsGraph based object to another.
  !> @{

  !> @brief Epetra_OffsetIndex generates and index of offsets allowing direct access to data for Import/Export operations on Epetra_CrsGraph based objects such as Epetra_CrsMatrix.

  type Epetra_OffsetIndex !,extends(universal)     :: Epetra_OffsetIndex !"shell"
  end type

contains
  !> @name Constructor Function
  !> @{

  !> @brief Constructs a Epetra_OffsetIndex object from the graphs and an exporter.
  type(Epetra_OffsetIndex) function Epetra_OffsetIndex(SourceGraph,TargetGraph,exporter)
   class(Epetra_CrsGraph),intent(in) :: SourceGraph
   class(Epetra_CrsGraph),intent(in) :: TargetGraph
   type(Epetra_Export), intent(in) :: exporter
  end function

  !> @name Constructor Function
  !> @{

  !> @brief Constructs a Epetra_OffsetIndex object from the graphs and an importer.
  type(Epetra_OffsetIndex) function Epetra_OffsetIndex(SourceGraph,TargetGraph,importer)
   class(Epetra_CrsGraph),intent(in) :: SourceGraph
   class(Epetra_CrsGraph),intent(in) :: TargetGraph
   type(Epetra_Import), intent(in) :: importer
  end function

  !> @name Constructor Function
  !> @{

   !> @brief Epetra_OffsetIndex copy constructor.
  type(Epetra_OffsetIndex) function Epetra_OffsetIndex(this)
    type(Epetra_OffsetIndex) ,intent(in) :: this
  end function

end module 
