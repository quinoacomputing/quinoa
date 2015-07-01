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

module FEpetra_Export
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_Export_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_Comm  ,only: Epetra_Comm
  use FEpetra_BlockMap ,only: Epetra_BlockMap
  use iso_c_binding ,only: c_int
  use forepetra
  implicit none
  !private                   ! Hide everything by default
  !public :: Epetra_Export ! Expose type/constructors/methods

   !> <BR> Epetra_Export: This class builds an export object for efficient exporting of off-processor elements.

   !> @brief Epetra_Export is used to construct a communication plan that can be called repeatedly by computational classes such the Epetra matrix, vector and multivector classes to efficiently send data to a target processor.
    !! This class currently has one constructor, taking two Epetra_Map or Epetra_BlockMap objects.  The first map specifies the global IDs that are owned by the calling processor.  The second map specifies the global IDs of elements that we want to export to later.
 
  type Epetra_Export !,extends(universal)   :: Epetra_Export !"shell"
  contains
     ! Public member functions
     !procedure        :: NumSameIDs
     !procedure        :: NumPermuteIDs
     !procedure        :: NumRemoteIDs
     !procedure        :: NumExportIDs
     !procedure        :: NumSend
     !procedure        :: NumRecv
     !procedure        :: SourceMap
     !procedure        :: TargetMap
  end type

contains
  ! Original C++ prototype:
  ! Epetra_Export( const Epetra_BlockMap & SourceMap, const Epetra_BlockMap & TargetMap );
  ! CTrilinos prototype:
  ! CT_Epetra_Export_ID_t Epetra_Export_Create ( CT_Epetra_BlockMap_ID_t SourceMapID, CT_Epetra_BlockMap_ID_t TargetMapID );
 
  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_Export Constructor
  !> @brief Constructs a Epetra_Export object from the source and target maps.
  !! This constructor builds an Epetra_Export object by comparing the GID lists of the source and target maps.
  type(Epetra_Export) function Epetra_Export(Source_Map,Target_Map)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t
    class(Epetra_BlockMap), intent(in) :: Source_Map &
   !< SourceMap (In) Map containing the GIDs from which data should be exported from each processor to the target map whenever an export operation is performed using this exporter.
    class(Epetra_BlockMap), intent(in) :: Target_Map &
   !< TargetMap (In) Map containing the GIDs that should be used for exporting data.
   !< \n
   !< Warning: Note that the TargetMap \e must have GIDs uniquely owned, each GID of the target map can occur only lds an export object that will transfer objects built with SourceMap to objects built with TargetMap.
   !< \n\n
   !<  A Epetra_Export object categorizes the elements of the target map into three sets as follows:
   !< \n
   !< \li All elements in the target map that have the same GID as the corresponding element of the source map, starting with the first element in the target map, going up to the first element that is different from the source map. The number of these IDs is returned by NumSameIDs().
   !< \li All elements that are local to the processor, but are not part of the first set of elements. These elements have GIDs that are owned by the calling processor, but at least the first element of this list is permuted. Even if subsequent elements are not permuted, they are included in this list.  The number of permuted elements is returned by NumPermutedIDs().  The list of elements (local IDs) in the source map that are permuted can be found in the list PermuteFromLIDs().  The list of elements (local IDs) in the target map that are the new locations of the source elements can be found in the list PermuteToLIDs().
    !< \li All remaining elements of the target map correspond to global IDs that are owned by remote processors.  The number of these elements is returned by NumRemoteIDs() and the list of these is returned by RemoteLIDs().
   !< \n See Trilinos Epetra_Import documentation for an example.
  end function

  ! Original C++ prototype:
  ! Epetra_Export(const Epetra_Export& Exporter);
  ! CTrilinos prototype:
  ! CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( CT_Epetra_Export_ID_t ExporterID );

  !> @name Constructor Functions  
  !! @{  

  !> <BR> Constructs a copy of a Epetra_Export object
  !> @brief This constructor returns a copy of an Epetra_Export object.
  type(Epetra_Export) function Epetra_Export(this)
    type(Epetra_Export) ,intent(in) :: this
  end function

  !> @name Public Member Functions
  !! @{  

  !> @brief
  !! Returns the number of elements that are identical between the source and target maps, up to the first different ID.
  integer(c_int) function NumSameIDs(this)
    class(Epetra_Export), intent(in) :: this
  end function

  !> @name Public Member Functions
  !! @{  

  !> @brief
  !! Returns the number of elements that are local to the calling processor, but not part of the first NumSameIDs() elements.
  integer(c_int) function NumPermuteIDs(this)
    class(Epetra_Export), intent(in) :: this
  end function

  !> @name Public Member Functions
  !! @{  

  !> @brief
  !! Returns the number of elements that are not on the calling processor.
  integer(c_int) function NumRemoteIDs(this)
    class(Epetra_Export), intent(in) :: this
  end function

  !> @name Public Member Functions
  !! @{  

  !> @brief
  !! Returns the number of elements that must be sent by the calling processor to other processors.
  integer(c_int) function NumExportIDs(this)
    class(Epetra_Export), intent(in) :: this
  end function

  !> @name Public Member Functions
  !! @{  

  !> @brief
  !! Total number of elements to be sent.
  integer(c_int) function NumSend(this)
    class(Epetra_Export), intent(in) :: this
  end function

  !> @name Public Member Functions
  !! @{  

  !> @brief
  !! Total number of elements to be received.
  integer(c_int) function NumRecv(this)
    class(Epetra_Export), intent(in) :: this
  end function

  !> @name Public Member Functions
  !! @{  

  !> @brief
  !! Returns the SourceMap used to construct this exporter.
  type(Epetra_BlockMap) function SourceMap(this)
   class(Epetra_Export), intent(in) :: this
  end function

  !> @name Public Member Functions
  !! @{  

  !> @brief
  !! Returns the TargetMap used to construct this exporter.
  type(Epetra_BlockMap) function TargetMap(this)
   class(Epetra_Export), intent(in) :: this
  end function
end module 
