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

module FEpetra_Map
  use ForTrilinos_enums ,only: FT_Epetra_BlockMap_ID,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_error
  use FEpetra_Comm       ,only: Epetra_Comm
  use FEpetra_BlockMap   ,only: Epetra_BlockMap
  use iso_c_binding      ,only: c_int
  use forepetra
  implicit none
  private                     ! Hide everything by default
  public :: Epetra_Map        ! Expose type/constructors/methods

  type , extends(Epetra_BlockMap) :: Epetra_Map 
  end type

contains
  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_Map constructor for a Epetra-defined uniform linear distribution of elements. 
  type(Epetra_Map) function Epetra_Map(Num_GlobalElements,IndexBase,comm)
    use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int)   ,intent(in) :: Num_GlobalElements &
      !< Number of elements to distribute.
    integer(c_int)   ,intent(in) :: IndexBase &
      !< Minimum index value used for arrays that use this map. Typically 1 for Fortran
    class(Epetra_Comm),intent(in) :: comm &
      !< Polymorphic class Epetra_Comm communicator containing information on the number of processors.
  end function

  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_Map constructor for a user-defined linear distribution of elements. 
  type(Epetra_Map) function Epetra_Map(Num_GlobalElements,Num_MyElements,IndexBase,comm)
    use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int)    ,intent(in) :: Num_GlobalElements &
     !< Number of elements to distribute.
    integer(c_int)    ,intent(in) :: Num_MyElements &
     !< Number of elements owned by the calling processor.
    integer(c_int)    ,intent(in) :: IndexBase  &
     !< Minimum index value used for arrays that use this map. Typically 1 for Fortran
    class(Epetra_Comm) ,intent(in) :: comm  &
     !< Polymorphic class Epetra_Comm communicator containing information on the number of processors.
  end function
  
  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_Map constructor for a user-defined arbitrary distribution of elements. 
  type(Epetra_Map) function Epetra_Map(Num_GlobalElements,Num_MyElements,My_GlobalElements,IndexBase,comm)
    use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in)              :: Num_GlobalElements &
     !< Number of elements to distribute. Must be either -1 or equal to the computed sum of NumMyElements across all processors in the Epetra_Comm communicator
    integer(c_int) ,intent(in)              :: Num_MyElements !< Number of elements owned by the calling processor.
    integer(c_int) ,intent(in) ,dimension(:),allocatable:: My_GlobalElements &
     !< Integer array of length NumMyElements. The ith entry contains the global index value of the ith element on this processor. Index values are not required to be contiguous on a processor, or to be within the range of 0 to NumGlobalElements. As long as the index values are consistently defined and used, any set of NumGlobalElements distinct integer values is acceptable.
    integer(c_int) ,intent(in)              :: IndexBase  &
     !< Minimum index value used for arrays that use this map. Typically 1 for Fortran
    class(Epetra_Comm)                      :: comm &
     !< Polymorphic class Epetra_Comm communicator containing information on the number of processors.
  end function

  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_Map copy constructor
  type(Epetra_Map) function Epetra_Map(this) 
    type(Epetra_Map) ,intent(in) :: this 
  end function

end module 

