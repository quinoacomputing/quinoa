! Runs case for 3 velocity components in 3D space
! u(x,y,z)=10.0sin(x)
! v(x,y,z)=0.0
! w(x,y,z)=0.0
program main
#include "ForTrilinos_config.h"
#ifdef HAVE_MPI
  use mpi
  use FEpetra_MpiComm   ,only: Epetra_MpiComm
#else
  use FEpetra_SerialComm,only: Epetra_SerialComm
#endif
  use ForTrilinos_utils ,only: valid_kind_parameters
  use iso_c_binding     ,only: c_int,c_double
  use field_module ,only: initial_field
  use vector_field_module,only: vector_field
  use initializer  ,only: u_3D_initial,v_3D_initial,w_3D_initial,zero
  implicit none
#ifdef HAVE_MPI
  type(Epetra_MpiComm) :: comm
#else
  type(Epetra_SerialComm) :: comm
#endif
  type(vector_field)                      :: u,N_u,u_half
  procedure(initial_field) ,pointer :: initial_u,initial_v,initial_w,initial
  real(c_double) :: dt,half=0.5,t=0.,t_final=0.45,nu=1.
  real(c_double) :: t_start,t_end
  integer(c_int) :: tstep=1
  integer(c_int) :: ierr 
 
  if (.not. valid_kind_parameters()) stop 'C interoperability not supported on this platform.'
#ifdef HAVE_MPI
  call MPI_INIT(ierr) 
  t_start=MPI_Wtime()
  comm = Epetra_MpiComm(MPI_COMM_WORLD)
#else
  call cpu_time(t_start) 
  comm = Epetra_SerialComm()
#endif
  initial_u => u_3D_initial
  initial_v => v_3D_initial
  initial_w => w_3D_initial
  initial => zero
  u = vector_field(initial_u,initial,initial,comm)
  N_u = vector_field(initial,initial,initial,comm)
  u_half = vector_field(initial,initial,initial,comm)
  do while (t<=t_final)  !2nd-order Runge-Kutta:
    dt = u%runge_kutta_stable_step(2,nu) 
    N_u = u.dotGradient.u
    u_half = u + ( (((.laplacian.u)*nu) - N_u)*(dt*half) ) ! first substep
    !u_half = u + ( (.laplacian.u)*nu - N_u )*(dt*half) ! first substep
    ! Additional paranthesis are needed for NAG to be able to construct 
    ! results for operators correctly
    N_u = u_half.dotGradient.u_half
    u  = u + ( (((.laplacian.u_half)*nu) - N_u)*dt ) ! second substep
    !u  = u + ( (.laplacian.u_half)*nu - N_u )*dt ! second substep
    ! Additional paranthesis are needed for NAG to be able to construct 
    ! results for operators correctly
    t = t + dt
    if (comm%MyPID()==0) print *,'timestep=',tstep,'time=',t
    tstep=tstep+1
  end do
 if (comm%MyPID()==0) write(10,*) 'u at t=',t
#ifdef HAVE_MPI
  t_end= MPI_Wtime()
#else
  call cpu_time(t_end)
#endif
  if (comm%MyPID()==0) write(10,*) 'Elapsed CPU time=',t_end-t_start
  call u%output()
  call N_u%force_finalize
  call u_half%force_finalize
  call u%force_finalize
  call comm%force_finalize
#ifdef HAVE_MPI
  call MPI_FINALIZE(ierr)
#endif
end program
