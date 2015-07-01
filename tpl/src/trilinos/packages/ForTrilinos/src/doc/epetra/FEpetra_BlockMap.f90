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


module FEpetra_BlockMap
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal ,only : universal
  use ForTrilinos_error ,only : error
  use FEpetra_Comm  ,only : Epetra_Comm
  use iso_c_binding ,only : c_int
  use forepetra
  implicit none
  !private                   ! Hide everything by default
  !public :: Epetra_BlockMap ! Expose type/constructors/methods

  type :: Epetra_BlockMap !,extends(universal)      :: Epetra_BlockMap 
  contains
     !I take out all type-bound procedure so there is no duppplicate in the     documnetation.  The second copy has a inheritance graph but doesn't list th    e arguments of the functions
    end type

 
contains

  !> @name Constructor Functions
  !! @{
 
  !> <BR> Epetra_BlockMap constructor for a Epetra-defined uniform linear distribution of constant size elements.
  !>  @brief Creates a map that distributes NumGlobalElements elements evenly across all processors in the Epetra_Comm communicator. If NumGlobalElements does not divide exactly into the number of processors, the first processors in the communicator get one extra element until the remainder is gone.
  !! The elements are defined to have a constant fixed size specified by ElementSize.
  type(Epetra_BlockMap) function Epetra_BlockMap(Num_GlobalElements,Element_Size,IndexBase,comm) 
    integer(c_int) ,intent(in) :: Num_GlobalElements &
     !< Number of elements to distribute.
    integer(c_int) ,intent(in) :: Element_Size &
     !< Number of points or vector entries per element.
    integer(c_int) ,intent(in) :: IndexBase &
     !< Minimum index value used for arrays that use this map. Typically 1 for Fortran.
    class(Epetra_Comm)         :: comm &
     !< Polymorphic type Epetra_Comm communicator containing information on the number of processors.
  end function

  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_BlockMap constructor for a user-defined linear distribution of constant size elements.
  !> @brief Creates a map that puts NumMyElements on the calling processor. If NumGlobalElements=-1, the number of global elements will be the computed sum of NumMyElements across all processors in the Epetra_Comm communicator.
  !!The elements are defined to have a constant fixed size specified by ElementSize.
  type(Epetra_BlockMap) function Epetra_BlockMap(Num_GlobalElements,Num_MyElements,Element_Size,IndexBase,comm)
    integer(c_int) ,intent(in) :: Num_GlobalElements &
     !< Number of elements to distribute. Must be either -1 or equal to the computed sum of NumMyElements across all processors in the Epetra_Comm communicator.
    integer(c_int) ,intent(in) :: Num_MyElements &
     !< Number of elements owned by the calling processor.
    integer(c_int) ,intent(in) :: Element_Size &
     !< Number of points or vector entries per element.
    integer(c_int) ,intent(in) :: IndexBase &
     !< Minimum index value used for arrays that use this map. Typically 1 for Fortran.
    class(Epetra_Comm)         :: comm &
     !< Polymorphic type Epetra_Comm communicator containing information on the number of processors.
  end function

  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_BlockMap constructor for a user-defined arbitrary distribution of constant size elements.
  !> @brief Creates a map that puts NumMyElements on the calling processor. The indices of the elements are determined from the list MyGlobalElements. If NumGlobalElements=-1, the number of global elements will be the computed sum of NumMyElements across all processors in the Epetra_Comm communicator.
  !!The elements are defined to have a constant fixed size specified by ElementSize.
  type(Epetra_BlockMap) function Epetra_BlockMap(Num_GlobalElements,&
                                                        My_GlobalElements,Element_Size,IndexBase,comm)
    integer(c_int) ,intent(in) :: Num_GlobalElements &
     !< Number of elements to distribute. Must be either -1 or equal to the computed sum of NumMyElements across all processors in the Epetra_Comm communicator.
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements &
     !< Integer array of length NumMyElements. The ith entry contains the global index value of the ith element on this processor. Index values are not required to be contiguous on a processor, or to be within the range of 1 to NumGlobalElements. As long as the index values are consistently defined and used, any set of NumGlobalElements distinct integer values is acceptable.
    integer(c_int) ,intent(in) :: Element_Size &
     !< Number of points or vector entries per element.
    integer(c_int) ,intent(in) :: IndexBase &
     !< Minimum index value used for arrays that use this map. Typically 1 for Fortran.
    class(Epetra_Comm)         :: comm &
     !< Polymorphic type Epetra_Comm communicator containing information on the number of processors.
  end function


  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_BlockMap constructor for a user-defined arbitrary distribution of variable size elements.
  !> @brief Creates a map that puts NumMyElements on the calling processor. If NumGlobalElements=-1, the number of global elements will be the computed sum of NumMyElements across all processors in the Epetra_Comm communicator. 
  !! The elements are defined to have a variable size defined by ElementSizeList.
  type(Epetra_BlockMap) function Epetra_BlockMap(Num_GlobalElements,&
                                                       My_GlobalElements,Element_SizeList,IndexBase,comm) 
    integer(c_int) ,intent(in) :: Num_GlobalElements &
     !< Number of elements to distribute. Must be either -1 or equal to the computed sum of NumMyElements across all processors in the Epetra_Comm communicator.
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements &
     !< Integer array of length NumMyElements. The ith entry contains the global index value of the ith element on this processor. Index values are not required to be contiguous on a processor, or to be within the range of 1 to NumGlobalElements. As long as the index values are consistently defined and used, any set of NumGlobalElements distinct integer values is acceptable.
    integer(c_int) ,intent(in) ,dimension(:) :: Element_SizeList & 
     !< A list of the element sizes for elements owned by the calling processor. The ith entry contains the element size of the ith element on this processor.  
    integer(c_int) ,intent(in) :: IndexBase  &
     !<  Minimum index value used for arrays that use this map. Typically 1 for Fortran.
    class(Epetra_Comm)         :: comm  &
     !< Polymorphic type Epetra_Comm communicator containing information on the number of processors.
  end function
  
  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_BlockMap copy constructor
  type(Epetra_BlockMap) function Epetra_BlockMap(this)
    type(Epetra_BlockMap) ,intent(in) :: this &
     !< Epetra_BlockMap object to copy 
  end function

  !> @name Size and dimension accessor functions
  !! @{

  !> <BR> Number of elements across all processors.
  integer(c_int) function NumGlobalElements(this) 
    type(Epetra_BlockMap) ,intent(in) :: this
  end function 
 
  !> @name Size and dimension accessor functions
  !! @{

  !> <BR> Number of elements on the calling processor. 
  integer(c_int) function NumMyElements(this)
    type(Epetra_BlockMap) ,intent(in) :: this
  end function 

  !> @name Size and dimension accessor functions
  !! @{

  !> <BR> Index base for this map. 
  integer(c_int) function IndexBase(this)
    type(Epetra_BlockMap) ,intent(in) :: this
  end function 

  !> @name Miscellaneous boolean tests
  !! @{

  !> <BR> Returns true if this and Map are identical maps.
  logical function  SameAs(lhs,rhs)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    type(Epetra_BlockMap)        ,intent(in) :: lhs
    type(Epetra_BlockMap)        ,intent(in) :: rhs
  end function SameAs

  !> @name Miscellaneous boolean tests
  !! @{

  !> <BR> Returns true if this and Map have identical point-wise structure.
  logical function  PointSameAs(lhs,rhs)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    type(Epetra_BlockMap)        ,intent(in) :: lhs
    type(Epetra_BlockMap)        ,intent(in) :: rhs
  end function PointSameAs

  !> @name Array accessor functions
  !! @{ 

  !> <BR> Array containing list of global IDs assigned to the calling processor. 
  function MyGlobalElements(this) result(MyGlobalElementsList)
    type(Epetra_BlockMap)     ,intent(in)    :: this
    integer(c_int),dimension(:),allocatable   :: MyGlobalElementsList
    integer(c_int)                            :: junk
  end function 

  !> @name Size and dimension accessor functions
  !! @{ 

  !> <BR> Returns the size of elements in the map; only valid if map has constant element size.
  integer(c_int) function ElementSize(this)
    type(Epetra_BlockMap) ,intent(in) :: this
  end function 

  !> @name Size and dimension accessor functions
  !! @{ 

  !> <BR> Size of element for specified L_ID
  integer(c_int) function ElementSize(this,L_ID)
    type(Epetra_BlockMap) ,intent(in) :: this
    integer(c_int)         ,intent(in) :: L_ID
  end function 

  !> @name Miscellaneous boolean tests
  !! @{

  !> <BR> Returns true if the global ID space is contiguously divided (but not necessarily uniformly) across all processors. 
  logical function LinearMap(this)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    type(Epetra_BlockMap),intent(in) :: this
  end function

  !> @name Miscellaneous boolean tests
  !! @{

  !> <BR> Returns true if map is defined across more than one processor. 
  logical function DistributedGlobal(this)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    type(Epetra_BlockMap),intent(in) :: this
  end function


end module 

