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


module FEpetra_CrsGraph
  use ForTrilinos_enums ,only: FT_Epetra_CrsGraph_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_enum_wrappers ,only: FT_Epetra_DataAccess_E_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_BlockMap  ,only: Epetra_BlockMap 
  use iso_c_binding     ,only: c_int,c_double,c_char
  use forepetra
  implicit none
  !private                      ! Hide everything by default
  !public :: Epetra_CrsGraph ! Expose type/constructors/methods

  !> <BR> Epetra_CrsGraph: A class for constructing and using sparse compressed row graphs.

  !> @brief Epetra_CrsGraph enables the piecewise construction and use of sparse matrix graphs (the integer structure without values) where entries are intended for row access.  
  !! Epetra_CrsGraph is an attribute of all Epetra row-based matrix classes, defining their nonzero structure and also holding their Epetra_Map attributes.

  type Epetra_CrsGraph  !,extends(universal)  :: Epetra_CrsGraph !"shell"
  contains
  end type
   
 contains
  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_CrsGraph constuctor 
  !> @brief Epetra_CrsGraph constructor with fixed number of indices per row
  !! Creates a Epetra_CrsGraph object and allocates storage.
  type(Epetra_CrsGraph) function Epetra_CsGraph(CV,RowMap,NumIndicesPerRow,StaticProfile)
   use ForTrilinos_enums ,only: FT_boolean_t,FT_TRUE,FT_FALSE
   use iso_c_binding     ,only: c_int
   integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV &
   !< (In) A FT_Epetra_DataAccess_E_t enumerated type set to Copy or View.
   class(Epetra_BlockMap) ,intent(in) :: RowMap &
   !< (In) An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap) listing the rows that this processor will contribute to.
   integer(c_int)         ,intent(in) :: NumIndicesPerRow &
   !< (In) An integer that indicates the (approximate if StaticProfile=false) number of entries in the each row. Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
   logical,        optional                   :: StaticProfile  &
   !< (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernes.
  end function

  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_CrsGraph constuctor 
  !> @brief Epetra_CrsGraph constructor with variable number of indices per row
  !! Creates a Epetra_CrsGraph object and allocates storage.
  type(Epetra_CrsGraph) function Epetra_CrsGraph(CV,RowMap,NumIndicesPerRow,StaticProfile)
   use ForTrilinos_enums ,only: FT_boolean_t,FT_TRUE,FT_FALSE
   use iso_c_binding     ,only: c_int
   integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV &
   !< (In) A FT_Epetra_DataAccess_E_t enumerated type set to Copy or View.
   class(Epetra_BlockMap) ,intent(in) :: RowMap &
   !< (In) An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap) listing the rows that this processor will contribute to
   integer(c_int),dimension(:) ,intent(in) :: NumIndicesPerRow &
   !< NumIndicesPerRow - (In) An integer array of length NumMyRows such that NumIndicesPerRow(i) indicates the (approximate if StaticProfile=false) number of entries in the ith row.
   logical,        optional                   :: StaticProfile  &
   !< StaticProfile - (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernels.
  end function

  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_CrsGraph copy constuctor 
  !> @brief This will create a Level 1 deep copy. 
  !!This Graph will share ownership of the CrsGraphData object with the right hand side Graph.
  type(Epetra_CrsGraph) function Epetra_CrsGraph(this)
    type(Epetra_CrsGraph) ,intent(in) :: this
  end function

end module 

