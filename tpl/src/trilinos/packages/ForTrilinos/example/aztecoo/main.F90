program main
#include "ForTrilinos_config.h"
#ifdef HAVE_MPI
  use mpi
  use FEpetra_MpiComm,only:Epetra_MpiComm
#else
  use FEpetra_SerialComm,only:Epetra_SerialComm
#endif
  use FEpetra_Map
  use FEpetra_Vector
  use FEpetra_CrsMatrix
  use FAztecOO
  use ForTrilinos_error
  use ForTrilinos_enum_wrappers
  use iso_c_binding, only : c_int,c_double
  implicit none
#ifdef HAVE_MPI
  type(Epetra_MpiComm) :: comm
#else
  type(Epetra_SerialComm) :: comm
#endif
  type(Epetra_Map) :: map
  type(Epetra_CrsMatrix) :: A
  type(Epetra_Vector) :: x,b
  type(AztecOO) :: Solver
  type(error) ::err
  integer(c_int), parameter :: grid_resolution=8
  integer :: rc,ierr 
  integer(c_int),parameter::IndexBase=1
  integer(c_int) :: NumGlobalElements = 8, NumMyElements, NumEntries, MaximumIter=100
  integer(c_int),dimension(:),allocatable :: MyGlobalElements, NumNonZero
  real(c_double),dimension(:),allocatable:: c
  real(c_double) :: values(2)
  real(c_double), parameter :: tolerance = 1.0E-4,three=3.0
  integer(c_int) :: MyGlobalElements_diagonal(1),i,indices(2)
  integer(c_int),parameter :: diagonal=1
  logical        :: success = .true.
 
#ifdef HAVE_MPI
  call MPI_INIT(ierr) 
  comm = Epetra_MpiComm(MPI_COMM_WORLD)
#else
  comm = Epetra_SerialComm()
#endif
  map = Epetra_Map(NumGlobalElements,IndexBase,comm)
  NumMyElements = map%NumMyElements()
  allocate(MyGlobalElements(NumMyElements))
  MyGlobalElements = map%MyGlobalElements()
  allocate(NumNonZero(NumMyElements))
  NumNonZero = 3  ! Non zero elements per row
! Create a Epetra_Matrix
  A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNonZero)
! Add rows one at a time
! off diagonal values will always be 1 and 1
  values(1) = 1.0
  values(2) = 1.0
  do i=1,NumMyElements
    if (MyGlobalElements(i)==1) then
      indices(1) = NumGlobalElements 
      indices(2) = 2
      NumEntries = 2
    else if(MyGlobalElements(i)==NumGlobalElements) then
      indices(1) = NumGlobalElements-1
      indices(2) = 1
      NumEntries = 2
    else
      indices(1) = MyGlobalElements(i)-1
      indices(2) = MyGlobalElements(i)+1
      NumEntries = 2
    end if
     call A%InsertGlobalValues(MyGlobalElements(i),values,indices)
  !Put in the diaogonal entry
     MyGlobalElements_diagonal=MyGlobalElements(i)
     call A%InsertGlobalValues(MyGlobalElements(i),[three],MyGlobalElements_diagonal)
  end do
  !Finish up
  call A%FillComplete(.true.)
  x=Epetra_Vector(map,.true.)
  b=Epetra_Vector(map,.true.)
  call x%Random()
  call b%PutScalar(1.0_c_double)
  Solver=AztecOO(A,x,b)
  call Solver%iterate(A,x,b,MaximumIter,tolerance,err)
  allocate(c(x%MyLength()))
  c=x%ExtractCopy(err)
  do i=1, NumMyElements
   if (abs((c(i)-0.2)/0.2)>tolerance) success = .false.
   print *,c(i)
  enddo 
   if (success) then
    print *
    print *, "End Result: TEST PASSED"
  else
    print *
    print *, "End Result: TEST FAILED"
  end if
  call Solver%force_finalize
  call A%force_finalize
  call b%force_finalize
  call x%force_finalize 
  call map%force_finalize
  call comm%force_finalize
#ifdef HAVE_MPI
  call MPI_FINALIZE(rc)
#endif
end program
