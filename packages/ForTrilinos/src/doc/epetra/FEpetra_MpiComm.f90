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
module FEpetra_MpiComm
  use ForTrilinos_enums ,only : FT_Epetra_Comm_ID,FT_Epetra_MpiComm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_error ,only : error
  use FEpetra_Comm      ,only : Epetra_Comm
  use iso_c_binding     ,only : c_int,c_long,c_double,c_char
  use forepetra
  implicit none
  !private                     ! Hide everything by default
  !public :: Epetra_MpiComm ! Expose type/constructors/methods

  type ,extends(Epetra_Comm)        :: Epetra_MpiComm 
  contains
    ! !Barrier Methods
    ! procedure         :: barrier
    ! !Broadcast Methods
    ! procedure :: broadcast_double
    ! procedure :: broadcast_int
    ! procedure :: broadcast_long
    ! procedure :: broadcast_char
    ! !Gather Methods
    ! procedure :: gather_double
    ! procedure :: gather_int
    ! procedure :: gather_long
    ! !Sum Methods
    ! procedure :: sum_double
    ! procedure :: sum_int
    ! procedure :: sum_long
    ! !Max/Min Methods
    ! procedure :: max_double
    ! procedure :: max_int
    ! procedure :: max_long
    ! procedure :: min_double
    ! procedure :: min_int
    ! procedure :: min_long
    ! !Parallel Prefix Methods
    ! procedure :: ScanSum_double
    ! procedure :: ScanSum_int
    ! procedure :: ScanSum_long
    ! !Attribute Accessor Methods
    ! procedure :: MyPID
    ! procedure :: NumProc
     !Gather/catter and Directory Constructors
     !I/O methods
  end type

contains
  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_MpiComm MPI Constructor.  
  !> @brief Creates a Epetra_MpiComm instance for use with MPI.  If no specialized MPI communicator is needed, this constuctor can be called with the argument MPI_COMM_WORLD.  
  type(Epetra_MpiComm) function Epetra_MpiComm(MPI_COMM_WORLD)
  end function
 
  !> @name Constructor Functions
  !! @{
  
  !> <BR> Epetra_MpiComm Copy Constructor. 
  !> @brief Makes an exact copy of an existing Epetra_MpiComm instance.
  type(Epetra_MpiComm) function Epetra_MpiComm(this)
    type(Epetra_MpiComm) ,intent(in) :: this 
  end function

  !> @name Barrier Methods 
  !! @{

  !> <BR> Epetra_MpiComm Barrier function. 
  !> @brief Causes each processor in the communicator to wait until all processors have arrived.
  !! Implements Epetra_Comm.
  subroutine barrier(this)
    class(Epetra_MpiComm) ,intent(in) :: this
  end subroutine
 
  !> @name Broadcast Methods
  !! @{

  !> <BR> Epetra_MpiComm Broadcast function. 
  !> @brief Takes list of input values from the root processor and sends to all other processors.
  !!  Implements Epetra_Comm.
  subroutine broadcast(this,MyVals,root,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(inout) :: MyVals &
     !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
    integer(c_int)               ,intent(in)    :: root &
     !< In On entry, contains the processor from which all processors will receive a copy of MyVals.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine
  
  !> @name Broadcast Methods
  !! @{

  !> <BR> Epetra_MpiComm Broadcast function. 
  !> @brief Takes list of input values from the root processor and sends to all other processors.
  !!  Implements Epetra_Comm.
  subroutine broadcast(this,MyVals,root,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(inout) :: MyVals &
     !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
    integer(c_int)               ,intent(in)    :: root &
     !< In On entry, contains the processor from which all processors will receive a copy of MyVals.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information
  end subroutine

  !> @name Broadcast Methods
  !! @{

  !> <BR> Epetra_MpiComm Broadcast function. 
  !> @brief Takes list of input values from the root processor and sends to all other processors.
  !!  Implements Epetra_Comm.
  subroutine broadcast_long(this,MyVals,root,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_long),dimension(:) ,intent(inout) :: MyVals &
     !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
    integer(c_int)               ,intent(in)    :: root &
     !< In On entry, contains the processor from which all processors will receive a copy of MyVals.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information
  end subroutine
 
  !> @name Broadcast Methods
  !! @{

  !> <BR> Epetra_MpiComm Broadcast function. 
  !> @brief Takes list of input values from the root processor and sends to all other processors.
  !!  Implements Epetra_Comm.
  subroutine broadcast(this,MyVals,root,err)
    class(Epetra_MpiComm)           ,intent(in)    :: this
    character(kind=c_char),dimension(:),intent(inout) :: MyVals &
     !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
    integer(c_int)                     ,intent(in)    :: root &
     !< In On entry, contains the processor from which all processors will receive a copy of MyVals.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine
  
  !> @name GatherAll Methods
  !! @{

  !> <BR> Epetra_MpiComm All Gather function. 
  !> @brief Take list of input values from all processors in the communicator and creates an ordered contiguous list of those values on each processor.
  !!  Implements Epetra_Comm.
 subroutine GatherAll(this,MyVals,AllVals,err)
   class(Epetra_MpiComm)     ,intent(in)    :: this
   real(c_double), dimension(:) ,intent(in)    :: MyVals &
   !< On entry, contains the list of values, to be sent to all processors.
   real(c_double), dimension(:) ,intent(inout) :: AllVals &
   !< On exit, contains the list of values from all processors. Must be of size NumProc*size(MyVals).
   type(error) ,optional, intent(inout) :: err &
   !< Return any error information.
  end subroutine

  !> @name GatherAll Methods
  !! @{

  !> <BR> Epetra_MpiComm All Gather function. 
  !> @brief Take list of input values from all processors in the communicator and creates an ordered contiguous list of tho    se values on each processor.
  !!  Implements Epetra_Comm.
  subroutine GatherAll(this,MyVals,AllVals,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: MyVals &
    !< On entry, contains the list of values, to be sent to all processors.
    integer(c_int), dimension(:) ,intent(inout) :: AllVals &
    !< On exit, contains the list of values from all processors. Must be of size NumProc*size(MyVals).
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  !> @name GatherAll Methods
  !! @{

  !> <BR> Epetra_MpiComm All Gather function. 
  !> @brief Take list of input values from all processors in the communicator and creates an ordered contiguous list of tho    se values on each processor.
  !!  Implements Epetra_Comm.
  subroutine gather_long(this,MyVals,AllVals,err)
    class(Epetra_MpiComm)      ,intent(in)    :: this
    integer(c_long), dimension(:) ,intent(in)    :: MyVals &
    !< On entry, contains the list of values, to be sent to all processors.
    integer(c_long), dimension(:) ,intent(inout) :: AllVals &
    !< On exit, contains the list of values from all processors. Must be of size NumProc*size(MyVals).
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  !> @name Sum Methods
  !! @{

  !> <BR> Epetra_MpiComm Global Summ function. 
  !> @brief Take list of input values from all processors in the communicator, computes the sum and returns the sum to all processors.
  !!  Implements Epetra_Comm.
  subroutine SumAll(this,PartialSums,GlobalSums,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialSums &
    !<  On entry, contains the list of values, usually partial sums computed locally, to be summed across all processors.
    real(c_double), dimension(:) ,intent(inout) :: GlobalSums &
    !<   On exit, contains the list of values summed across all processors.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.   
  end subroutine

  !> @name Sum Methods
  !! @{

  !> <BR> Epetra_MpiComm Global Summ function. 
  !> @brief Take list of input values from all processors in the communicator, computes the sum and returns the sum to all     processors.
  !!  Implements Epetra_Comm.
  subroutine SumAll(this,PartialSums,GlobalSums,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialSums &
    !<  On entry, contains the list of values, usually partial sums computed locally, to be summed across all processors.
    integer(c_int), dimension(:) ,intent(inout) :: GlobalSums &
    !<   On exit, contains the list of values summed across all processors.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine

  !> @name Sum Methods
  !! @{

  !> <BR> Epetra_MpiComm Global Summ function. 
  !> @brief Take list of input values from all processors in the communicator, computes the sum and returns the sum to all     processors.
  !!  Implements Epetra_Comm.
  subroutine sum_long(this,PartialSums,GlobalSums,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialSums &
    !<  On entry, contains the list of values, usually partial sums computed locally, to be summed across all processors.
    integer(c_long), dimension(:),intent(inout) :: GlobalSums &
    !<   On exit, contains the list of values summed across all processors.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine
 
  !> @name Max/Min Methods
  !! @{

  !> <BR> Epetra_MpiComm Global Max function. 
  !> @brief Take list of input values from all processors in the communicator, computes the max and returns the max to all processors.
  !!  Implements Epetra_Comm. 
  subroutine MaxAll(this,PartialMaxs,GlobalMaxs,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialMaxs &
    !<  On entry, contains the list of values, usually partial maxs computed locally, using these Partial Maxs, the max across all processors will be computed.
    real(c_double), dimension(:) ,intent(inout) :: GlobalMaxs &
    !<  On exit, contains the list of maxs computed across all processors.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  !> @name Max/Min Methods
  !! @{

  !> <BR> Epetra_MpiComm Global Max function. 
  !> @brief Take list of input values from all processors in the communicator, computes the max and returns the max to all     processors.
  !!  Implements Epetra_Comm. 
  subroutine MaxAll(this,PartialMaxs,GlobalMaxs,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialMaxs &
    !<  On entry, contains the list of values, usually partial maxs computed locally, using these Partial Maxs, the max across all processors will be computed.
    integer(c_int), dimension(:) ,intent(inout) :: GlobalMaxs &
    !< On exit, contains the list of maxs computed across all processors.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  !> @name Max/Min Methods
  !! @{

  !> <BR> Epetra_MpiComm Global Max function. 
  !> @brief Take list of input values from all processors in the communicator, computes the max and returns the max to all     processors.
  !!  Implements Epetra_Comm. 
  subroutine max_long(this,PartialMaxs,GlobalMaxs,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialMaxs &
    !<  On entry, contains the list of values, usually partial maxs computed locally, using these Partial Maxs, the max across all processors will be computed.
    integer(c_long), dimension(:),intent(inout) :: GlobalMaxs &
    !< On exit, contains the list of maxs computed across all processors.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine
 
  !> @name Max/Min Methods
  !! @{
    
  !> <BR> Epetra_MpiComm Global Min function. 
  !> @brief  Take list of input values from all processors in the communicator, computes the min and returns the min to all processors.
  !!  Implements Epetra_Comm. 
  subroutine MinAll(this,PartialMins,GlobalMins,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialMins &
    !<  On entry, contains the list of values, usually partial mins computed locally; using these Partial Mins, the min across all processors will be computed.
    real(c_double), dimension(:) ,intent(inout) :: GlobalMins &
    !< On exit, contains the list of mins computed across all processors.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  !> @name Max/Min Methods
  !! @{
    
  !> <BR> Epetra_MpiComm Global Min function. 
  !> @brief  Take list of input values from all processors in the communicator, computes the min and returns the min to all     processors.  
  !!  Implements Epetra_Comm. 
  subroutine MinAll(this,PartialMins,GlobalMins,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialMins &
    !<  On entry, contains the list of values, usually partial mins computed locally; using these Partial Mins, the min across all processors will be computed.
    integer(c_int), dimension(:) ,intent(inout) :: GlobalMins &
    !< On exit, contains the list of mins computed across all processors.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine

  !> @name Max/Min Methods
  !! @{

  !> <BR> Epetra_MpiComm Global Min function. 
  !> @brief  Take list of input values from all processors in the communicator, computes the min and returns the min to all     processors.
  !!  Implements Epetra_Comm.
  subroutine min_long(this,PartialMins,GlobalMins,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialMins &
    !<  On entry, contains the list of values, usually partial mins computed locally; using these Partial Mins, the min across all processors will be computed.
    integer(c_long), dimension(:),intent(inout) :: GlobalMins &
    !< On exit, contains the list of mins computed across all processors.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine

  !> @name Parallel Prefix Methods
  !! @{

  !> <BR> Epetra_MpiComm Scan Sum function. 
  !> @brief Take list of input values from all processors in the communicator, computes the scan sum and returns it to all processors such that processor i contains the sum of values from processor 0 up to and including processor i.
  !!  Implements Epetra_Comm.
  subroutine ScanSum(this,MyVals,scan_sums,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: MyVals  &
    !< On entry, contains the list of values to be summed across all processors.
    real(c_double), dimension(:) ,intent(inout) :: scan_sums &
    !< On exit, contains the list of values summed across processors 0 through i.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  !> @name Parallel Prefix Methods
  !! @{

  !> <BR> Epetra_MpiComm Scan Sum function. 
  !> @brief Take list of input values from all processors in the communicator, computes the scan sum and returns it to all     processors such that processor i contains the sum of values from processor 0 up to and including processor i.
  !!  Implements Epetra_Comm.
  subroutine ScanSum(this,MyVals,scan_sums,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: MyVals  &
    !< On entry, contains the list of values to be summed across all processors.
    integer(c_int), dimension(:) ,intent(inout) :: scan_sums &
    !< On exit, contains the list of values summed across processors 0 through i.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine

  !> @name Parallel Prefix Methods
  !! @{
    
  !> <BR> Epetra_MpiComm Scan Sum function.  
  !> @brief Take list of input values from all processors in the communicator, computes the scan sum and returns it to all     processors such that processor i contains the sum of values from processor 0 up to and including processor i.
  !!  Implements Epetra_Comm.
  subroutine ScanSum_long(this,MyVals,scan_sums,err)
    class(Epetra_MpiComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: MyVals &
    !< On entry, contains the list of values to be summed across all processors.
    integer(c_long), dimension(:),intent(inout) :: scan_sums &
    !< On exit, contains the list of values summed across processors 0 through i.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine
 
  !< name Attribute Accessor Methods
  !! @{
    
  !> <BR> Return my process ID.  
  !> @brief In MPI mode returns the rank of the calling process.  In serial mode returns 0.
  !!  Implements Epetra_Comm.
  integer(c_int) function MyPID(this)
    class(Epetra_MpiComm)     , intent(in) :: this
  end function
  
  !< name Attribute Accessor Methods
  !! @{

  !> <BR> Return my process ID.  
  !> @brief Returns total number of processes (always returns 1 for MpiComm).
  !!  Implements Epetra_Comm.
  integer(c_int) function NumProc(this)
    class(Epetra_MpiComm)     , intent(in) :: this
  end function

end module 
