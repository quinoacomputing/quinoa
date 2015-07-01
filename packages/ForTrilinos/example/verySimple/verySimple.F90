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
  use iso_c_binding ,only : c_int,c_double,c_bool
  use iso_fortran_env ,only : error_unit ,output_unit
  use fortrilinos_utils ,only : valid_kind_parameters
  use forepetra 
  implicit none

  ! This file is the Fortran equivalent of CTrilionos/example/verySimple.c.
  ! As such, it makes direct use of the procedural bindings in forepetra.F90.
  ! Although we deprecate the use of this style, we retain this file for instructive
  ! purposes: users who are learning Fortran 2003 and know an earlier version of Fortran
  ! might find it useful to compare this Fortran 90-style code with the analogous object-
  ! oriented code in ../verySimpleObjectOriented.F90.
  
  !
  ! Data declarations 
  !

  integer(c_int) numGlobalElements, numGlobalElements_rtn ;
  integer(c_int) junk;

  type(FT_Epetra_Comm_ID_t) commID;

  type(FT_Epetra_BlockMap_ID_t) bmapID;

  type(FT_Epetra_Vector_ID_t) xID, bID;
  type(FT_Epetra_MultiVector_ID_t) mxID, mbID;

  real(c_double) bnorm(1), xnorm(1), expected_bnorm, expected_xnorm, bnorm_err, xnorm_err ;
  real(c_double) ,parameter :: err_tol=0.0

  logical(c_bool) :: success = .TRUE.;

  !
  ! Executable code
  !

  if (.not. valid_kind_parameters()) &
     stop 'In ForTrilinos (verySimple.F90): C interoperability not supported on this platform.'
  
  ! Create an Epetra_SerialComm but store it as an Epetra_Comm so that
  ! it can be passed to functions expecting the latter 

  commID = Epetra_Comm_Degeneralize(Epetra_SerialComm_Generalize(Epetra_SerialComm_Create()));

  ! Create an Epetra_Map and store it as an Epetra_BlockMap so that
  ! a) it can be passed to functions expecting the latter and
  ! b) methods implemented only in BlockMap can be invoked on the Map

  numGlobalElements = 4;
  ! use indexBase = 0 unless you know what you're doing!
  bmapID = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(Epetra_Map_Create(numGlobalElements, indexBase=0,CommID=commID)));

  ! Check the properties of the map 
  numGlobalElements_rtn = Epetra_BlockMap_NumGlobalElements(bmapID);
  print *, "NumGlobalElements = ", numGlobalElements_rtn ;
  if ( numGlobalElements /= numGlobalElements_rtn ) stop 'main: incorrect numGlobalElements_rtn';

  print *,"NumMyElements = ", Epetra_BlockMap_NumMyElements(bmapID);
  
  ! Create an Epetra_Vector and store it as an Epetra_MultiVector so that
  ! methods implemented only in MultiVector can be invoked on the Vector 
  ! zero this one 
  mxID = Epetra_MultiVector_Degeneralize(Epetra_Vector_Generalize(Epetra_Vector_Create(bmapID, FT_TRUE)));

  ! Do the same thing, but do not initialize this one to zero 
  mbID = Epetra_MultiVector_Degeneralize(Epetra_Vector_Generalize(Epetra_Vector_Create(bmapID, FT_FALSE)));

  ! Do some vector operations 
  junk = Epetra_MultiVector_PutScalar(mbID, 2.0_c_double);
  junk = Epetra_MultiVector_Update_WithA(mxID, 2.0_c_double, mbID, 0.0_c_double);!/* x = 2*b */

  junk = Epetra_MultiVector_Norm2(mbID, bnorm);
  junk = Epetra_MultiVector_Norm2(mxID, xnorm);

  print *,"2 norm of x = ", xnorm;
  print *,"2 norm of b = ", bnorm;

  ! Test the expected value 
  expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements );
  expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements );
  bnorm_err = abs( expected_bnorm - bnorm(1) ) / expected_bnorm;
  xnorm_err = abs( expected_xnorm - xnorm(1) ) / expected_xnorm;
  print *,"error in 2 norm of x = ", bnorm_err ;
  print *,"error in 2 norm of b = ", xnorm_err ;

  if (bnorm_err > err_tol) success = .FALSE.;
  if (xnorm_err > err_tol) success = .FALSE.;

  ! Clean up memory (in reverse order)
  call Epetra_MultiVector_Destroy(mxID);
  call Epetra_MultiVector_Destroy(mbID);
  call Epetra_BlockMap_Destroy(bmapID);
  call Epetra_Comm_Destroy(commID);

  if (success) then
    write(output_unit,*) 
    write(output_unit,fmt='(a)') "End Result: TEST PASSED" ;
  else
    write(error_unit,*) 
    write(error_unit,fmt='(a)') "End Result: TEST FAILED" ;
  end if
  
end program
