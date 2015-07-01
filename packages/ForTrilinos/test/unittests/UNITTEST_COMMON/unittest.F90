program main
#include "ForTrilinos_config.h"

  use iso_fortran_env ,only : error_unit ,output_unit
  use ForTrilinos_external_utils
  use TEST_CALLS_FILE
#ifdef HAVE_MPI
  use mpi
#endif
  implicit none

  logical :: success,fullsuccess,multi
  character(len=100) :: which_test
  character(len=100),dimension(100) :: test_list
  integer :: ierr,test_num,test_cnt,rank

#ifdef HAVE_MPI
  call MPI_INIT(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
#else
  rank = 0
#endif

  success = .FALSE.
  call get_test_list(rank, test_list, test_cnt, multi)
#ifdef HAVE_MPI
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

  fullsuccess = .TRUE.
  do test_num = 1,test_cnt
    which_test = test_list(test_num)

    if (rank == 0) print *,"Testing ",TEST_IMPLS_FILE_STR, &
        & "::",trim(which_test),"_UnitTest"

    call test_setup()
    success = select_test(which_test)
    call test_teardown()

#ifdef HAVE_MPI
    call reduce_success(success)
#endif

    if (rank == 0) then
      if (success) then
        write(output_unit,fmt='(a)') "TEST PASSED"
      else
        write(output_unit,fmt='(a)') "TEST FAILED"
      end if
    end if

    fullsuccess = fullsuccess .and. success
    flush(output_unit)
  end do

  if (multi) then
    if (rank == 0) then
      if (fullsuccess) then
        write(output_unit,fmt='(a)') "END RESULT: ALL TESTS PASSED" ;
      else
        write(output_unit,fmt='(a)') "END RESULT: SOME TESTS FAILED" ;
      end if
      flush(output_unit)
    end if
  end if

#ifdef HAVE_MPI
  call MPI_FINALIZE(ierr)
#endif

  contains

    subroutine read_test_list(fname, test_list, test_cnt)
      character(len=100),intent(in) :: fname
      character(len=100),dimension(100),intent(out) :: test_list
      integer,intent(out) :: test_cnt

      integer :: ierr,u,test_num
      character(len=100) :: which_test

      u = 100
      open(unit=u, iostat=ierr, status='old', file=fname)
      if (ierr .NE. 0) call handle_err("Error opening test list file "//fname, __LINE__)
      test_cnt = 0
      do test_num = 1,100
        read(u, '(a100)', iostat=ierr) which_test
        if (ierr .NE. 0) call handle_err("Error in test list file "//which_test, __LINE__)
        if (which_test == "###ENDOFTESTS") exit
        test_list(test_num) = which_test
        test_cnt = test_cnt + 1
      end do
      close(unit=u)
    end subroutine

#ifndef HAVE_MPI
    subroutine get_test_list(rank,list,cnt,multi)
      integer,intent(in) :: rank
      character(len=100),dimension(100),intent(out) :: list
      integer,intent(out) :: cnt
      logical,intent(out) :: multi

      character(len=100) :: which_test
      integer :: test_num,arg_len,stat,arg_cnt

      arg_cnt = command_argument_count()
      if (arg_cnt == 0) call handle_err("No tests specified", __LINE__)

      call get_command_argument(1,which_test,arg_len,stat)
      if (stat < 0) call handle_err("Command line argument cropped", __LINE__)
      if (stat > 0) call handle_err("Could not retrieve first command line argument", __LINE__)

      if (which_test == "-f") then
        if (arg_cnt < 2) call handle_err("No test list file specified", __LINE__)
        call get_command_argument(2,which_test,arg_len,stat)
        if (stat < 0) call handle_err("Test file command line argument cropped", __LINE__)
        if (stat > 0) call handle_err("Could not retrieve command line argument for test file", __LINE__)
        call read_test_list(which_test, list, cnt)
        multi = .TRUE.
      else
        cnt = arg_cnt
        if (cnt > 100) call handle_err("Too many tests specified", __LINE__)
        do test_num = 1,cnt
          call get_command_argument(test_num,which_test,arg_len,stat)
          test_list(test_num) = which_test
        end do
        if (cnt > 1) then
          multi = .TRUE.
        else
          multi = .FALSE.
        end if
      end if
    end subroutine
#else
    subroutine get_test_list(rank,list,cnt,multi)
      integer,intent(inout) :: rank
      character(len=100),dimension(100),intent(out) :: list
      integer,intent(out) :: cnt
      logical,intent(out) :: multi

      integer,dimension(1) :: passcnt
      character(len=100) :: arg
      integer :: test_num,arg_len,stat,arg_cnt,ierr

      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

      if (rank == 0) then
        arg_cnt = command_argument_count()
        if (arg_cnt == 0) call handle_err("No tests specified", __LINE__)
        passcnt(1) = arg_cnt
      end if
      call MPI_Bcast(passcnt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      arg_cnt = passcnt(1)

      if (rank == 0) then
        call get_command_argument(1,arg,arg_len,stat)
        if (stat < 0) call handle_err("Command line argument cropped", __LINE__)
        if (stat > 0) call handle_err("Could not retrieve first command line argument", __LINE__)
      end if
      call MPI_Bcast(arg, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      if (arg == "-f") then
        if (arg_cnt < 2) call handle_err("No test list file specified", __LINE__)
        if (rank == 0) then
          call get_command_argument(2,arg,arg_len,stat)
          if (stat < 0) call handle_err("Test file command line argument cropped", __LINE__)
          if (stat > 0) call handle_err("Could not retrieve command line argument for test file", __LINE__)
        end if
        call MPI_Bcast(arg, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call read_test_list(arg, list, cnt)
        multi = .TRUE.
      else
        call handle_err("You must specify a test file in MPI mode",__LINE__)
      end if
    end subroutine
#endif

#ifdef HAVE_MPI
    subroutine reduce_success(success)
      logical,intent(inout) :: success

      integer :: ierr
      integer,dimension(1) :: localfail,anyfail

      if (success) then
        localfail(1) = 0
      else
        localfail(1) = 1
      end if
      call MPI_AllReduce(localfail, anyfail, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (anyfail(1) == 0) then
        success = .TRUE.
      else
        success = .FALSE.
      end if
    end subroutine
#endif

    subroutine test_setup()
#ifdef HAVE_MPI
      print *, "Starting with fresh tables..."
      call ForTrilinos_CleanSlate()
#endif
    end subroutine

    subroutine test_teardown()
    end subroutine

    subroutine handle_err(msg, line)
      character(len=*),intent(in) :: msg
      integer,intent(in) :: line

      integer :: ierr

      print *, "ERROR: "//trim(msg)//" (line ", line, &
              & " of "//__FILE__//"). TEST FAILED"
      flush(output_unit)
#ifdef HAVE_MPI
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
      stop 1
#endif
    end subroutine
  
end program
