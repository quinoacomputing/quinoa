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
module FEpetra_Comm
  use ForTrilinos_universal ,only : universal
  use ForTrilinos_enums !,only: FT_Epetra_Comm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_error
  use ForTrilinos_table_man
  use forepetra
  implicit none
  !private               ! Hide everything by default
  !public :: Epetra_Comm ! Expose type/methods

  !> <BR> Epetra_Comm:  The Epetra Communication Abstract Base Class.
  !> @{

  !> @brief The Epetra_Comm class is an interface that encapsulates the general information and services needed for other Epetra classes to run on a parallel computer. An Epetra_Comm object is required for building all Epetra Map objects, which in turn are required for all other Epetra classes.
  !! Epetra_Comm has default implementations, via Epetra_SerialComm and Epetra_MpiComm, for both serial execution and MPI distributed memory execution.  It is meant to insulate the user from the specifics of communication that are not required for normal manipulation of linear algebra objects.  

  type, abstract :: Epetra_Comm !,abstract ,extends(universal) :: Epetra_Comm
  contains
    !Barrier Methods
  !  procedure(barrier)          ,deferred  ::barrier
    !Broadcast Methods
  !  procedure(broadcast) ,deferred  ::broadcast_double
  !  procedure(broadcast)    ,deferred  ::broadcast_int
  !  procedure(broadcast_long_interface)   ,deferred  ::broadcast_long
  !  procedure(broadcast_char_interface)   ,deferred  ::broadcast_char
  !  generic :: broadcast=>broadcast_double,broadcast_int,broadcast_char
  !  !Gather Methods
  !  procedure(gather_double_interface),deferred  ::gather_double
  !  procedure(gather_int_interface)   ,deferred  ::gather_int
  !  procedure(gather_long_interface)  ,deferred  ::gather_long
  !  generic :: GatherAll=>gather_double,gather_int
  !  !Sum Methods
  !  procedure(sum_double_interface),deferred   ::sum_double
  !  procedure(sum_int_interface)   ,deferred   ::sum_int
  !  procedure(sum_long_interface)  ,deferred   ::sum_long
  !  generic :: SumAll=>sum_double,sum_int
  !  !Max/Min Methods
  !  procedure(max_double_interface) ,deferred   ::max_double
  !  procedure(max_int_interface)    ,deferred   ::max_int
  !  procedure(max_long_interface)   ,deferred   ::max_long
  !  generic :: MaxAll=>max_double,max_int
  !  procedure(min_double_interface) ,deferred   ::min_double
  !  procedure(min_int_interface)    ,deferred   ::min_int
  !  procedure(min_long_interface)               ,deferred   ::min_long
  !  generic :: MinAll=>min_double,min_int
  !  !Parallel Prefix Methods
  !  procedure(ScanSum_double_interface) ,deferred   ::ScanSum_double
  !  procedure(ScanSum_int_interface)    ,deferred   ::ScanSum_int
  !  procedure(ScanSum_long_interface)   ,deferred   ::ScanSum_long
  !  generic :: ScanSum=>ScanSum_double,ScanSum_int
  !  !Attribute Accessor Methods
  !  procedure(MyPID_interface)           ,deferred::MyPID
  !  procedure(NumProc_interface)         ,deferred::NumProc
  !  !Gather/catter and Directory Constructors
  !  !I/O methods
  end type
  
  abstract interface

    !> @name Barrier Methods
    !! @{
 
    !> <BR> Epetra_Comm Barrier function. 
    !> @brief Each processor must wait at the point the barrier is called until all processors have arrived.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine barrier(this) 
      import:: Epetra_Comm
      class(Epetra_Comm) ,intent(in)  :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information 
    end subroutine

    !> @name Broadcast Methods
    !! @{

    !> <BR> Epetra_Comm Broadcast function.
    !> @brief Take list of input values from the root processor and sends to all other processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine broadcast(this,MyVals,root,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information 
      real(c_double) ,dimension(:) ,intent(inout) :: MyVals &
       !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
      integer(c_int)               ,intent(in)    :: root &
       !< In On entry, contains the processor from which all processors will receive a copy of Values.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Broadcast Methods
    !! @{

    !> <BR> Epetra_Comm Broadcast function.
    !> @brief Take list of input values from the root processor and sends to all other processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine broadcast(this,MyVals,root,err) 
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm,error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_int) ,dimension(:) ,intent(inout) :: MyVals &
        !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
      integer(c_int)               ,intent(in)    :: root &
       !< In On entry, contains the processor from which all processors will receive a copy of Values.
      type(error)   ,optional      ,intent(inout) :: err &
      !< Returns  error information.
    end subroutine
    
    !> @name Broadcast Methods
    !! @{

    !> <BR> Epetra_Comm Broadcast function.
    !> @brief Take list of input values from the root processor and sends to all other processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine broadcast_long(this,MyVals,root,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm,error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_long),dimension(:) ,intent(inout) :: MyVals &
       !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
      integer(c_int)               ,intent(in)    :: root&
       !< In On entry, contains the processor from which all processors will receive a copy of Values.
      type(error)   ,optional      ,intent(inout) :: err &
      !< Returns  error information.
    end subroutine
    
    !> @name Broadcast Methods
    !! @{

    !> <BR> Epetra_Comm Broadcast function.
    !> @brief Take list of input values from the root processor and sends to all other processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine broadcast(this,MyVals,root,err) 
      use iso_c_binding ,only: c_int,c_char
      import:: Epetra_Comm, error
      class(Epetra_Comm)                 ,intent(in)    :: this &
        !< Polymorphic type Epetra_Comm communicator containing processors information
      character(kind=c_char),dimension(:),intent(inout) :: MyVals &
        !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
      integer(c_int)                     ,intent(in)    :: root &
       !< In On entry, contains the processor from which all processors will receive a copy of Values.
      type(error)   ,optional            ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine
    
    !> @name Gather Methods
    !! @{

    !> <BR> Epetra_Comm All Gather function.
    !> @brief Take list of input values from all processors in the communicator and creates an ordered contiguous list of those values on each processor.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine GatherAll(this,MyVals,AllVals,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)         ,intent(in)    :: this &
        !< Polymorphic type Epetra_Comm communicator containing processors information
      real(c_double),dimension(:),intent(in)    :: MyVals &
        !< In On entry, contains the list of values to be sent to all processors.
      real(c_double),dimension(:),intent(inout) :: AllVals  &
        !< Out On exit, contains the list of values from all processors. Must be of size NumProc*size(MyVals).
      type(error)   ,optional    ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Gather Methods
    !! @{

    !> <BR> Epetra_Comm All Gather function.
    !> @brief Take list of input values from all processors in the communicator and creates an ordered contiguous list of those values on each processor.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine GatherAll(this,MyVals,AllVals,err) 
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)         ,intent(in)   :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_int),dimension(:),intent(in)   :: MyVals &
       !< In On entry, contains the list of values to be sent to all processors.
      integer(c_int),dimension(:),intent(inout):: AllVals &
       !< Out On exit, contains the list of values from all processors. Must be of size NumProc*size(MyVals).
      type(error)   ,optional    ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Gather Methods
    !! @{

    !> <BR> Epetra_Comm All Gather function.
    !> @brief Take list of input values from all processors in the communicator and creates an ordered contiguous list of those values on each processor.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine gather_long(this,MyVals,AllVals,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_long),dimension(:) ,intent(in)    :: MyVals &
       !< In On entry, contains the list of values to be sent to all processors.
      integer(c_long),dimension(:) ,intent(inout):: AllVals &
       !< Out On exit, contains the list of values from all processors. Must be of size NumProc*size(MyVals).
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Sum Methods
    !! @{

    !> <BR> Epetra_Comm Sum function.
    !> @brief Take list of input values from all processors in the communicator, computes the sum and returns the sum to all processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine SumAll(this,PartialSums,GlobalSums,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      real(c_double),dimension(:)  ,intent(in)    :: PartialSums &
       !< In On entry, contains the list of values, usually partial sums computed locally, to be summed across all processors.
      real(c_double),dimension(:)  ,intent(inout) :: GlobalSums &
       !< Out On exit, contains the list of values summed across all processors.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Sum Methods
    !! @{

    !> <BR> Epetra_Comm Sum function.
    !> @brief Take list of input values from all processors in the communicator, computes the sum and returns the sum to all processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine SumAll(this,PartialSums,GlobalSums,err) 
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
        !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_int), dimension(:) ,intent(in)    :: PartialSums &
        !< In On entry, contains the list of values, usually partial sums computed locally, to be summed across all processors.
      integer(c_int), dimension(:) ,intent(inout) :: GlobalSums &
       !< Out On exit, contains the list of values summed across all processors.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Sum Methods
    !! @{

    !> <BR> Epetra_Comm Sum function.
    !> @brief Take list of input values from all processors in the communicator, computes the sum and returns the sum to all processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine sum_long(this,PartialSums,GlobalSums,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_long), dimension(:),intent(in)    :: PartialSums &
       !< In On entry, contains the list of values, usually partial sums computed locally, to be summed across all processors.
      integer(c_long), dimension(:),intent(inout) :: GlobalSums &
       !< Out On exit, contains the list of values summed across all processors.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Max/Min Methods
    !! @{

    !> <BR> Epetra_Comm Global Max function.
    !> @brief Take list of input values from all processors in the communicator, computes the max and returns the max to all processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine MaxAll(this,PartialMaxs,GlobalMaxs,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      real(c_double), dimension(:) ,intent(in)    :: PartialMaxs &
       !< In On entry, contains the list of values, usually partial maxs computed locally; using these Partial Maxs, the max across all processors will be computed.
      real(c_double), dimension(:) ,intent(inout) :: GlobalMaxs &
       !< Out On exit, contains the list of maxs computed across all processors.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Max/Min Methods
    !! @{

    !> <BR> Epetra_Comm Global Max function.
    !> @brief Take list of input values from all processors in the communicator, computes the max and returns the max to all processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine MaxAll(this,PartialMaxs,GlobalMaxs,err) 
      use iso_c_binding ,only: c_int,c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_int), dimension(:) ,intent(in)    :: PartialMaxs &
       !< In On entry, contains the list of values, usually partial maxs computed locally; using these Partial Maxs, the max across all processors will be computed.
      integer(c_int), dimension(:) ,intent(inout) :: GlobalMaxs &
       !< Out On exit, contains the list of maxs computed across all processors.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine
    !> @name Max/Min Methods
    !! @{

    !> <BR> Epetra_Comm Global Max function.
    !> @brief Take list of input values from all processors in the communicator, computes the max and returns the max to all processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine max_long(this,PartialMaxs,GlobalMaxs,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_long), dimension(:),intent(in)    :: PartialMaxs &
       !< In On entry, contains the list of values, usually partial maxs computed locally; using these Partial Maxs, the max across all processors will be computed.
      integer(c_long), dimension(:),intent(inout) :: GlobalMaxs &
       !< Out On exit, contains the list of maxs computed across all processors.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Max/Min Methods
    !! @{

    !> <BR> Epetra_Comm Global Min function.
    !> @brief Take list of input values from all processors in the communicator, computes the min and returns the min to all processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine MinAll(this,PartialMins,GlobalMins,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      real(c_double), dimension(:) ,intent(in)    :: PartialMins &
       !< In On entry, contains the list of values, usually partial mins computed locally; using these Partial Mins, the min across all processors will be computed.
      real(c_double), dimension(:) ,intent(inout) :: GlobalMins &
      !< Out On exit, contains the list of mins computed across all processors.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Max/Min Methods
    !! @{

    !> <BR> Epetra_Comm Global Min function.
    !> @brief Take list of input values from all processors in the communicator, computes the min and returns the min to all processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine MinAll(this,PartialMins,GlobalMins,err) 
      use iso_c_binding ,only: c_int,c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
      !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_int), dimension(:) ,intent(in)    :: PartialMins &
      !< In On entry, contains the list of values, usually partial mins computed locally; using these Partial Mins, the min across all processors will be computed.
      integer(c_int), dimension(:) ,intent(inout) :: GlobalMins &
      !< Out On exit, contains the list of mins computed across all processors.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Max/Min Methods
    !! @{

    !> <BR> Epetra_Comm Global Min function.
    !> @brief Take list of input values from all processors in the communicator, computes the min and returns the min to all processors.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine min_long(this,PartialMins,GlobalMins,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
      !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_long), dimension(:),intent(in)    :: PartialMins &
      !< In On entry, contains the list of values, usually partial mins computed locally; using these Partial Mins, the min across all processors will be computed.
      integer(c_long), dimension(:),intent(inout) :: GlobalMins &
      !< Out On exit, contains the list of mins computed across all processors.
      type(error)   ,optional      ,intent(inout) :: err &
       !< Returns  error information.
    end subroutine

    !> @name Parallel Prefix Methods
    !! @{

    !> <BR> Epetra_Comm Scan Sum function.
    !> @brief Take list of input values from all processors in the communicator, computes the scan sum and returns it to all processors such that processor i contains the sum of values from processor 0 up to and including processor i.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine ScanSum(this,MyVals,scan_sums,err) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      real(c_double), dimension(:) ,intent(in)    :: MyVals &
       !< In On entry, contains the list of values to be summed across all processors.
      real(c_double), dimension(:) ,intent(inout) :: scan_sums &
       !< Out On exit, contains the list of values summed across processors 0 through i. 
      type(error)   ,optional      ,intent(inout) :: err  &
       !< Returns  error information.
    end subroutine

    !> @name Parallel Prefix Methods
    !! @{

    !> <BR> Epetra_Comm Scan Sum function.
    !> @brief Take list of input values from all processors in the communicator, computes the scan sum and returns it to all processors such that processor i contains the sum of values from processor 0 up to and including processor i.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine ScanSum(this,MyVals,scan_sums,err) 
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
      !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_int), dimension(:) ,intent(in)    :: MyVals &
      !< In On entry, contains the list of values to be summed across all processors.
      integer(c_int), dimension(:) ,intent(inout) :: scan_sums &
      !< Out On exit, contains the list of values summed across processors 0 through i.
      type(error)   ,optional      ,intent(inout) :: err &
      !< Returns  error information.
    end subroutine

    !> @name Parallel Prefix Methods
    !! @{

    !> <BR> Epetra_Comm Scan Sum function.
    !> @brief Take list of input values from all processors in the communicator, computes the scan sum and returns it to all processors such that processor i contains the sum of values from processor 0 up to and including processor i.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    subroutine ScanSum_long(this,MyVals,scan_sums,err) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm, error
      class(Epetra_Comm)           ,intent(in)    :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_long), dimension(:),intent(in)    :: MyVals &
      !< In On entry, contains the list of values to be summed across all processors.
      integer(c_long), dimension(:),intent(inout) :: scan_sums &
      !< Out On exit, contains the list of values summed across processors 0 through i.
      type(error)   ,optional      ,intent(inout) :: err &
      !< Returns  error information.
    end subroutine

    !> @name Attribute Accessor Methods
    !! @{

    !> <BR> Return my process ID.
    !> @brief In MPI mode returns the rank of the calling process. In serial mode returns 0.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    function MyPID(this)
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm
      class(Epetra_Comm), intent(in) :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_int) :: MyPID
    end function

    !> @name Attribute Accessor Methods
    !! @{

    !> <BR> Returns total number of processes.
    !> @brief In MPI mode returns the size of the MPI communicator. In serial mode returns 1.
    !!Implemented in Epetra_MpiComm, and Epetra_SerialComm.
    function NumProc(this)
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm
      class(Epetra_Comm), intent(in) :: this &
       !< Polymorphic type Epetra_Comm communicator containing processors information
      integer(c_int) :: NumProc
    end function
  end interface


end module 
