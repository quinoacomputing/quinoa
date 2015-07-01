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

program main
  ! This file is the object-oriented fortran equivalent of Epetra_power_method.cpp. In Trilinos 10.4,
  ! this is a snapshot of an unstable (evolving) file expected to become stable in a
  ! subsequent release.  This file exercises the derived types defined in 
  ! ForTrilinos/src/epetra/FEpetra*.F90, which wrap the interface bodies in 
  ! ForTrilinos/src/epetra/forepetra.F90.   
    
  ! This file represents the preferred style for using ForTrilinos and is recommended for 
  ! Fortran users whose compilers support the object-oriented features of Fortran 2003.
  ! As of the Trilinos 10.4 release date, the latest versions of the IBM and Cray compilers 
  ! nominally support the required features.  The Numerical Algorithms Group (NAG) 
  ! supports all of the required features of Fortran 2003.  

#include "ForTrilinos_config.h"
#ifdef HAVE_MPI
  use mpi
  use FEpetra_MpiComm      ,only : Epetra_MpiComm
#else
  use FEpetra_SerialComm   ,only : Epetra_SerialComm
#endif
  use FEpetra_Map          ,only : Epetra_Map
  use FEpetra_Vector       ,only : Epetra_Vector
  use FEpetra_CrsMatrix    ,only : Epetra_CrsMatrix
  use ForTrilinos_utils    ,only : valid_kind_parameters
  use ForTrilinos_enum_wrappers
  use ForTrilinos_error
  use iso_c_binding        ,only : c_int,c_double
  implicit none
! Data declarations 
#ifdef HAVE_MPI  
  type(Epetra_MpiComm) :: communicator
#else
  type(Epetra_SerialComm) :: communicator
#endif
  type(Epetra_Map)    :: map
  type(Epetra_CrsMatrix) :: A
  type(error)         :: err
  integer(c_int)      :: NumGlobalElements, NumGlobalElements_return
  integer(c_int),dimension(:),allocatable :: MyGlobalElements
  integer(c_int),dimension(:),allocatable :: NumNz
  integer(c_int)      :: NumMyElements,i
  integer(c_int)      :: Index_Base=1
  integer(c_int)      :: MyPID, NumProc
  logical             :: verbose
  integer(c_int)      :: indices(2)
  real(c_double)      :: two = 2.0,values(2)
  integer             :: rc, ierr,ierr_pm 
  real(c_double)      :: lambda(1),tolerance=1.0E-2
  integer(c_int)      :: niters, numvals,index_diagonal
  integer(c_int),dimension(:),allocatable::Rowinds
  real(c_double),dimension(:),allocatable::Rowvals
  logical        :: success = .true.

  if (.not. valid_kind_parameters()) stop 'C interoperability not supported on this platform.'
  
  ! Executable code
  
! Create a comm
#ifdef HAVE_MPI
  call MPI_INIT(ierr)
  communicator= Epetra_MpiComm(MPI_COMM_WORLD) 
# else
  communicator= Epetra_SerialComm() 
#endif

  MyPID   = communicator%MyPID()
  NumProc = communicator%NumProc()
  verbose = MyPID==0
  print *,verbose
  print *,MyPID,NumProc

! Create a map 
  NumGlobalElements = 100 
  if (NumGlobalElements < NumProc) stop 'Number of global elements cannot be less that number of processors'
  map = Epetra_Map(NumGlobalElements,Index_Base,communicator)
  NumGlobalElements_return = map%NumGlobalElements()   ! test line

! Get updated list and number of local equations from newly created Map
  NumMyElements = map%NumMyElements()
  print *,'NumGlobalElements = ', numGlobalElements_return  ! test line
  print *,'NumMyElements=', map%NumMyElements()             ! test line
  if ( NumGlobalElements /= NumGlobalElements_return ) &
    stop 'In ForTrilinos (verySimpleObjectOriented.F90: return mismatch'
  allocate(MyGlobalElements(NumMyElements))
  MyGlobalElements = map%MyGlobalElements()

! Create an integer vector NumNz that is used to build the Epetra Matrix
! NumNz(i) is the number of OFF-DIAGONAL term for the ith global equation
! on this processor
  allocate(NumNz(NumMyElements))

! We are building a tridiagonal matrix where each row has (-1 2 -1)
! So we need 2 off-diagonal terms (except for the first and last equation)
  do i=1,NumMyElements
   if(MyGlobalElements(i)==1.or.MyGlobalElements(i)==NumGlobalElements) then
     NumNz(i) = 2
   else
     NumNz(i) = 3
   end if
  end do
  
! Create a Epetra_Matrix
  A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

! Add rows one at a time
! Need some vectors to help
! off diagonal values will always be -1
  values(1) = -1.0
  values(2) = -1.0
  do i=1,NumMyElements
    if (MyGlobalElements(i)==1) then
      indices(1) = 2
      call A%InsertGlobalValues(MyGlobalElements(i),[values(1)],[indices(1)],err)
    else if(MyGlobalElements(i)==NumGlobalElements) then
      indices(1) = NumGlobalElements-1
      call A%InsertGlobalValues(MyGlobalElements(i),[values(1)],[indices(1)],err)
    else
      indices(1) = MyGlobalElements(i)-1
      indices(2) = MyGlobalElements(i)+1
      call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
    end if
     if (err%error_code()/=0) stop 'A%InsertGlobalValues: failed'
  !Put in the diaogonal entry
     index_diagonal = MyGlobalElements(i)
     call A%InsertGlobalValues(MyGlobalElements(i),[two],[index_diagonal],err)
     if (err%error_code()/=0) stop 'A%InsertGlobalValues: failed'
  end do
 
  !Finish up
  call A%FillComplete()
 
  !Create vectors for power methods
  !variable needed for interation
  lambda = 0.0
  ierr_pm = 0
  niters=NumGlobalElements*10 
  
  !Iterate
  call power_method(A,lambda,niters,tolerance,verbose,ierr_pm)
  ierr_pm=ierr_pm+1

  if (A%MyGlobalRow(1_c_int)) then
    numvals=A%NumGlobalEntries(1_c_int)
    call A%ExtractGlobalRowCopy(1_c_int,Rowvals,Rowinds) ! get A(1,1)
  do i=1,numvals
   if (Rowinds(i)==1) Rowvals(i)=10.0*Rowvals(i)
  enddo  
   call A%ReplaceGlobalValues(1_c_int,Rowvals,Rowinds)
  endif

  !Iterate again
  lambda = 0.0
  call power_method(A,lambda,niters,tolerance,verbose,ierr_pm)
  ierr_pm=ierr_pm+1

 
  ! Clean up memory (in reverse order).  This step is not required
  ! with compilers that support Fortran 2003 type finalization:
  call A%force_finalize()
  call map%force_finalize()
  call communicator%force_finalize()
 
#ifdef HAVE_MPI
  call MPI_FINALIZE(rc)
#endif
  ! Should really figure out cases where this does not work
  if (success) then
    print *  
    print *, "End Result: TEST PASSED" 
  else
    print *  
    print *, "End Result: TEST FAILED"
  end if

contains

subroutine power_method(A,lambda,niters,tolerance,verbose,ierr_pm)
 implicit none
 type(Epetra_CrsMatrix), intent(inout) :: A 
 integer(c_int),    intent(inout):: ierr_pm
 integer(c_int), intent(in) :: niters
 real(c_double), intent(in) :: tolerance 
 logical, intent(in)        :: verbose
 real(c_double)             :: lambda(1)
 real(c_double) :: normz(1),residual(1)
 type(Epetra_Vector) :: q,z,resid
 integer(c_int) :: iter
 ierr_pm=1
 q = Epetra_Vector(A%RowMap()) 
 z = Epetra_Vector(A%RowMap()) 
 resid = Epetra_Vector(A%RowMap()) 

 !Fill z with random numbers
 call z%Random()

 do iter=0,niters
  normz=z%Norm2()              !Compute 2-norm of z
  call q%Scale(1.0/normz(1),z)  
  call A%Multiply(.false.,q,z) ! Compute z=A*q
  lambda = q%Dot(z)            ! Approximate maximum eignvalue
  if (modulo(iter,100)==0.or.iter+1==niters) then
   call resid%Update(1.0_c_double,z,-lambda(1),q,0.0_c_double) ! Compute A*q-lambda*q
   residual=resid%Norm2()
   if (verbose) print *,'Iter=',iter,'lambda=',lambda(1),'Resisual of A*q-lambda*q=',residual(1)
  endif
  if (residual(1)<tolerance) then
   ierr_pm=0
   exit
  endif
 enddo
end subroutine

end program main
