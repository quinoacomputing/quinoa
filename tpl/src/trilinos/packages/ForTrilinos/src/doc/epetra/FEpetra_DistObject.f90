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
module FEpetra_DistObject
  use ForTrilinos_enums ,only : FT_Epetra_SrcDistObject_ID,FT_Epetra_DistObject_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
  use ForTrilinos_universal
  use ForTrilinos_table_man
  use ForTrilinos_error ,only : error
  use FEpetra_SrcDistObject ,only : Epetra_SrcDistObject
  use FEpetra_Export, only: Epetra_Export
  use FEpetra_Import, only: Epetra_Import
  use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
  use iso_c_binding     ,only : c_int,c_long,c_double,c_char
  use forepetra
  implicit none
  !private                     ! Hide everything by default
  !public :: Epetra_DistObject ! Expose type/constructors/methods

  !> <BR> Epetra_DistObject: A class for constructing and using dense multi-vectors, vectors and matrices in parallel.
  !> @{

  !> @brief The Epetra_DistObject is a base class for all Epetra distributed global objects.  
  !! It provides the basic mechanisms and interface specifications for importing and exporting operations using Epetra_Import and Epetra_Export objects.

  type ,extends(Epetra_SrcDistObject)        :: Epetra_DistObject !"shell"
  contains
     !Import/Export methods
     !procedure, private:: DistObject_Export
     !procedure, private:: DistObject_Export_UsingImporter
     !generic :: export => DistObject_Export_UsingImporter, DistObject_Export
     !procedure, private:: DistObject_Import
     !procedure, private:: DistObject_Import_UsingExporter
     !generic :: import => DistObject_Import_UsingExporter, DistObject_Import
  end type

contains

  !> @name Import/Export Methods
  !> @{
 
  !> @brief Exports an Epetra_DistObject using the Epetra_Export object.
  subroutine export(this,A,exporter,CombineMode,indexor,err)
   class(Epetra_DistObject), intent(in) :: this 
   class(Epetra_SrcDistObject), intent(in) :: A &
   !< Distributed object that will be exported to the \e this multivector.
   type(Epetra_Export),intent(in) :: exporter &
   !< A Epetra_Export object specifying the communication required.
   integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode &
   !< A FT_Epetra_CombineMode_E_t enumerated type specifying how results should be combined on the receiving processor.
   type(Epetra_OffsetIndex), intent(in) :: indexor
   type(error),optional,intent(out) :: err
   !< Returns  error information.
  end subroutine

  !> @name Import/Export Methods
  !> @{
 
  !> @brief Exports an Epetra_DistObject using the Epetra_Import object.
  subroutine export(this,A,importer,CombineMode,indexor,err)
   class(Epetra_DistObject), intent(in) :: this
   class(Epetra_SrcDistObject), intent(in) :: A &
   !< Distributed object that will be exported to the \e this object.
   type(Epetra_Import),intent(in) :: importer &
   !< A Epetra_Import object specifying the communication required.
   integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode &
   !< A FT_Epetra_CombineMode_E_t enumerated type specifying how results should be combined on the receiving processor.
   type(Epetra_OffsetIndex), intent(in) :: indexor
   type(error),optional,intent(out) :: err &
   !< Returns  error information.
  end subroutine

  !> @name Import/Export Methods
  !> @{
 
  !> @brief Imports an Epetra_DistObject using the Epetra_Import object.
  subroutine import(this,A,importer,CombineMode,indexor,err)
   class(Epetra_DistObject), intent(in) :: this
   class(Epetra_SrcDistObject), intent(in) :: A &
   !< Distributed object that will be imported into the \e this object.
   type(Epetra_Import),intent(in) :: importer &
   !< A Epetra_Import object specifying the communication required.
   integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode &
   !< A FT_Epetra_CombineMode_E_t enumerated type specifying how results should be combined on the receiving processor.
   type(Epetra_OffsetIndex), intent(in) :: indexor
   type(error),optional,intent(out) :: err &
   !< Returns  error information.
  end subroutine

  !> @name Import/Export Methods
  !> @{
 
  !> @brief Imports an Epetra_DistObject using the Epetra_Export object.
  subroutine import(this,A,exporter,CombineMode,indexor,err)
   class(Epetra_DistObject), intent(in) :: this
   class(Epetra_SrcDistObject), intent(in) :: A &
   !< Distributed object that will be imported into the \e this object.
   type(Epetra_Export),intent(in) :: exporter &
   !< A Epetra_Export object specifying the communication required.
   integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode &
   !< A FT_Epetra_CombineMode_E_t enumerated type specifying how results should be combined on the receiving processor.
   type(Epetra_OffsetIndex), intent(in) :: indexor
   type(error),optional,intent(out) :: err &
   !< Returns  error information.
  end subroutine
end module 
