module initializer
  use iso_c_binding, only: c_double
  implicit none
  abstract interface
    real(c_double) pure function initial_field(x)
      import :: c_double
      real(c_double) ,intent(in) :: x
    end function
  end interface
contains
  real(c_double) pure function u_initial(x)
    real(c_double) ,intent(in) :: x
    u_initial = 10._c_double*sin(x)
  end function
  real(c_double) pure function zero(x)
    real(c_double) ,intent(in) :: x
    zero = 0.
  end function
end module 
module periodic_2nd_order_module
  use FEpetra_Comm, only:Epetra_Comm
  use FEpetra_Map!, only: Epetra_Map
  use FEpetra_Vector!, only:Epetra_Vector
  use FEpetra_CrsMatrix, only:Epetra_CrsMatrix
  use ForTrilinos_enum_wrappers
  use ForTrilinos_error 
  use initializer ,only : initial_field
  use iso_c_binding ,only : c_double,c_int
  use ForTrilinos_assertion_utility ,only: error_message,assert,assert_identical
  implicit none
  private
  public :: periodic_2nd_order
  type   :: periodic_2nd_order
    private
    type(Epetra_Vector) :: f
  contains
    procedure :: add => total
    procedure :: subtract => difference
    procedure :: multiply_field => product
    procedure :: multiply_real => multiple
    procedure :: runge_kutta_2nd_step => rk2_dt
    procedure :: x  => df_dx     ! 1st derivative w.r.t. x
    procedure :: xx => d2f_dx2   ! 2nd derivative w.r.t. x
    generic :: operator(+)   => add 
    generic :: operator(-)   => subtract
    generic :: operator(*)   => multiply_real,multiply_field
    procedure :: output
    procedure :: force_finalize
  end type
   
  real(c_double) ,parameter                 :: pi=acos(-1._c_double)
  real(c_double) ,dimension(:) ,allocatable :: x_node
  type(Epetra_Map)             ,allocatable :: map

  interface periodic_2nd_order
    procedure constructor
  end interface 

contains
  function constructor(initial,grid_resolution,comm) result(this)
    type(periodic_2nd_order) ,allocatable :: this
    procedure(initial_field) ,pointer :: initial
    integer(c_int) ,intent(in) :: grid_resolution
    integer(c_int) :: i
    class(Epetra_Comm), intent(in) :: comm
    integer(c_int) :: NumGlobalElements
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int)      :: NumMyElements,IndexBases=1,status
    real(c_double) ,dimension(:) ,allocatable :: f_v
    type(error) :: ierr
     
    NumGlobalElements=grid_resolution
    allocate(this)
    if (.not. allocated(x_node)) x_node = grid()
    if (.not. allocated(map)) then
      allocate(map,stat=status)
      ierr=error(status,'periodic_2nd_order: create map')
      call ierr%check_success()
      map = Epetra_Map(NumGlobalElements,IndexBases,comm)
    end if

    NumMyElements= map%NumMyElements()
    allocate(MyGlobalElements(NumMyElements))
    MyGlobalElements = map%MyGlobalElements()
    allocate(f_v(NumMyElements))
    forall(i=1:NumMyElements) f_v(i)=initial(x_node(MyGlobalElements(i))) 
    this%f=Epetra_Vector(map,zero_initial=.true.)
    call this%f%ReplaceGlobalValues(f_v,MyGlobalElements)
  contains
    pure function grid()
      integer(c_int) :: i
      real(c_double) ,dimension(:) ,allocatable :: grid
      allocate(grid(grid_resolution))
      forall(i=1:grid_resolution) &
        grid(i)  = 2.*pi*real(i-1,c_double)/real(grid_resolution,c_double)  
    end function
  end function

  real(c_double) function rk2_dt(this,nu, grid_resolution)
    class(periodic_2nd_order) ,intent(in) :: this
    real(c_double) ,intent(in) :: nu
    integer(c_int) ,intent(in) :: grid_resolution
    real(c_double)             :: dx, CFL, k_max
    dx=2.0*pi/grid_resolution
    k_max=grid_resolution/2.0_c_double
    CFL=1.0/(1.0-cos(k_max*dx))
    rk2_dt = CFL*dx**2/nu
  end function

  function total(lhs,rhs)
    class(periodic_2nd_order) ,intent(in) :: lhs
    class(periodic_2nd_order) ,intent(in) :: rhs
    type(periodic_2nd_order)  :: total
    total%f=Epetra_Vector(map,zero_initial=.true.)
    call total%f%Update(1._c_double,lhs%f,1._c_double,rhs%f,0._c_double)
  end function

  function difference(lhs,rhs)
    class(periodic_2nd_order) ,intent(in) :: lhs
    class(periodic_2nd_order) ,intent(in) :: rhs
    type(periodic_2nd_order) :: difference
    difference%f=Epetra_Vector(map,zero_initial=.true.)
    call difference%f%Update(1._c_double,lhs%f,-1._c_double,rhs%f,0._c_double)
  end function

 function product(lhs,rhs)
   class(periodic_2nd_order) ,intent(in) :: lhs
   class(periodic_2nd_order) ,intent(in) :: rhs
   type(periodic_2nd_order)  :: product
   product%f=Epetra_Vector(map,zero_initial=.true.)
   call product%f%Multiply(1._c_double,lhs%f,rhs%f,0._c_double)
end function

  function multiple(lhs,rhs)
    class(periodic_2nd_order) ,intent(in) :: lhs
    real(c_double) ,intent(in)  :: rhs
    type(periodic_2nd_order) :: multiple
    multiple%f=Epetra_Vector(map,zero_initial=.true.)
    call multiple%f%Scale(rhs,lhs%f)
  end function

  function df_dx(this) 
    class(periodic_2nd_order) ,intent(in) :: this
    type(periodic_2nd_order) :: df_dx
    type(Epetra_Vector) :: x
    type(Epetra_CrsMatrix) :: A 
    type(error) :: err
    real(c_double) :: dx
    integer(c_int) :: nx
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int) :: MyGlobalElements_diagonal(1)
    integer(c_int),dimension(:),allocatable :: NumNz
    integer(c_int) :: NumGlobalElements
    integer(c_int) :: NumMyElements,i
    integer(c_int) :: indices(2)
    real(c_double) ::values(2)
    real(c_double),parameter :: zero =0.0

  ! Executable code
   nx=size(x_node)
   dx=2.*pi/real(nx,c_double)
   NumGlobalElements = nx

! Get update list and number of local equations from given Map
  NumMyElements = map%NumMyElements()
  call assert_identical( [NumGlobalElements,map%NumGlobalElements()] )
  allocate(MyGlobalElements(NumMyElements))
  MyGlobalElements = map%MyGlobalElements()

! Create an integer vector NumNz that is used to build the Epetra Matrix
! NumNz(i) is the number of non-zero elements for the ith global equation
! on this processor
  allocate(NumNz(NumMyElements))

! We are building a tridiagonal matrix where each row has (-1 0 1)
! So we need 2 off-diagonal terms (except for the first and last equation)
  NumNz = 3

! Create a Epetra_Matrix
  A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

! Add rows one at a time
! Need some vectors to help
! off diagonal values will always be -1 and 1
  values(1) = -1.0/(2.0*dx)
  values(2) = 1.0/(2.0*dx)
  do i=1,NumMyElements
    if (MyGlobalElements(i)==1) then
      indices(1) = NumGlobalElements 
      indices(2) = 2
    else if(MyGlobalElements(i)==NumGlobalElements) then
      indices(1) = NumGlobalElements-1
      indices(2) = 1
    else
      indices(1) = MyGlobalElements(i)-1
      indices(2) = MyGlobalElements(i)+1
    end if
     call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
     call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
  !Put in the diaogonal entry
     MyGlobalElements_diagonal=MyGlobalElements(i)
     call A%InsertGlobalValues(MyGlobalElements(i),[zero],MyGlobalElements_diagonal,err)
     call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
  end do

  !Finish up
    call A%FillComplete(.true.,err)
    call assert( [err%error_code()==0_c_int] , [error_message('A%FillComplete: failed')] )
  !create vector x 
    x=Epetra_Vector(A%RowMap())
    Call A%Multiply_Vector(.false.,this%f,x)
 !create vector of df_dx
    df_dx%f=Epetra_Vector(x)
  end function
  
  function d2f_dx2(this) 
    class(periodic_2nd_order) ,intent(in) :: this
    type(periodic_2nd_order)  :: d2f_dx2
    type(Epetra_Vector) :: x
    type(Epetra_CrsMatrix) :: A 
    type(error) :: err
    real(c_double) :: dx
    integer(c_int) :: nx
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int) :: MyGlobalElements_diagonal(1)
    integer(c_int),dimension(:),allocatable :: NumNz
    integer(c_int) :: NumGlobalElements
    integer(c_int) :: NumMyElements,i
    integer(c_int) :: indices(2)
    real(c_double) :: values(2)
    real(c_double) :: two_dx2  

  ! Executable code
   nx=size(x_node)
   dx=2.*pi/real(nx,c_double)
   NumGlobalElements = nx

! Get update list and number of local equations from given Map
  NumMyElements = map%NumMyElements()
  call assert_identical( [NumGlobalElements,map%NumGlobalElements()] )
  allocate(MyGlobalElements(NumMyElements))
  MyGlobalElements = map%MyGlobalElements()

! Create an integer vector NumNz that is used to build the Epetra Matrix
! NumNz(i) is the number of non-zero elements for the ith global equation
! on this processor
  allocate(NumNz(NumMyElements))

! We are building a tridiagonal matrix where each row has (1 -2  1)
! So we need 2 off-diagonal terms (except for the first and last equation)
  NumNz = 3

! Create a Epetra_Matrix
  A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

! Add rows one at a time
! Need some vectors to help
! off diagonal values will always be 1 and 1
  values(1) = 1.0/(dx*dx)
  values(2) = 1.0/(dx*dx)
  two_dx2   =-2.0/(dx*dx)
  do i=1,NumMyElements
    if (MyGlobalElements(i)==1) then
      indices(1) = NumGlobalElements 
      indices(2) = 2
    else if(MyGlobalElements(i)==NumGlobalElements) then
      indices(1) = NumGlobalElements-1
      indices(2) = 1
    else
      indices(1) = MyGlobalElements(i)-1
      indices(2) = MyGlobalElements(i)+1
    end if
     call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
     call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
  !Put in the diaogonal entry
     MyGlobalElements_diagonal=MyGlobalElements(i)
     call A%InsertGlobalValues(MyGlobalElements(i),[two_dx2],MyGlobalElements_diagonal,err)
     call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
  end do

  !Finish up
    call A%FillComplete(.true.,err)
    call assert( [err%error_code()==0_c_int] , [error_message('A%FillComplete: failed')] )
  !create vector x 
    x=Epetra_Vector(A%RowMap())
    Call A%Multiply_Vector(.false.,this%f,x)
  !create vector of df_dx
    d2f_dx2%f=Epetra_Vector(x)
  end function

  subroutine output(this,comm)
    class(periodic_2nd_order) ,intent(in) :: this
    class(Epetra_Comm),intent(in) ::comm
    integer(c_int) :: i,NumMyElements,NumGlobalElements
    integer(c_int), dimension(:), allocatable :: MyGlobalElements
    real(c_double), dimension(:), allocatable :: f_v
    real(c_double), dimension(:), allocatable :: f
    NumGlobalElements=map%NumGlobalElements()
    NumMyElements=map%NumMyElements()
    allocate(MyGlobalElements(NumMyElements))
    MyGlobalElements=map%MyGlobalElements()
    allocate(f_v(NumMyElements))
    f_v=this%f%ExtractCopy()
    allocate(f(NumGlobalElements))
    call comm%GatherAll(f_v,f)
    do i=1,NumGlobalElements
     if (comm%MyPID()==0) write(20,'(2(E20.12,1x))') x_node(i),f(i)
    enddo
    call check_solution(f,NumGlobalElements)
  end subroutine

  subroutine check_solution(f,NumGlobalElements)
    integer(c_int) :: NumGlobalElements
    real(c_double), dimension(:), allocatable :: f,exact
    real(c_double),parameter :: tolerance=10.0E-1
    integer(c_int),parameter :: exact_grid=16
    logical :: success=.true.
    allocate(exact(exact_grid))
    exact=(/0.0, 0.686, 1.37, 2.053, 2.73, 3.387, 3.918, 3.57, & !exact solution at t=0.4626377
            0.0,-3.57,-3.918,-3.387,-2.73,-2.053,-1.37,-0.686/) !exact solution at t=0.4626377
    if (NumGlobalElements==exact_grid) then
     if  (sqrt(sum((f-exact)**2))>tolerance) success=.false.
    endif
    if (success) then
      print *
      print *, "End Result: TEST PASSED"
    else
      print *
      print *, "End Result: TEST FAILED"
    end if
  end subroutine 
  
  subroutine force_finalize(this)
    class(periodic_2nd_order), intent(inout) :: this
    call this%f%force_finalize
  end subroutine
end module
program main
#include "ForTrilinos_config.h"
#ifdef HAVE_MPI
  use mpi
  use FEpetra_MpiComm,only:Epetra_MpiComm
#else
  use FEpetra_SerialComm,only:Epetra_SerialComm
#endif
  use iso_c_binding, only : c_int,c_double
  use periodic_2nd_order_module, only : periodic_2nd_order
  use initializer ,only : u_initial,zero,initial_field
  implicit none
#ifdef HAVE_MPI
  type(Epetra_MpiComm) :: comm
#else
  type(Epetra_SerialComm) :: comm
#endif
  type(periodic_2nd_order), allocatable :: u,half_uu,u_half
  procedure(initial_field) ,pointer :: initial 
  real(c_double) :: dt,half=0.5,t=0.,t_final=0.4,nu=1.
  real(c_double) :: t_start,t_end
  integer(c_int) :: tstep=1
  integer(c_int), parameter :: grid_resolution=16
  integer :: rc,ierr 
 
#ifdef HAVE_MPI
  call MPI_INIT(ierr) 
  t_start=MPI_Wtime()
  comm = Epetra_MpiComm(MPI_COMM_WORLD)
#else
  call cpu_time(t_start) 
  comm = Epetra_SerialComm()
#endif
  initial => u_initial
  u = periodic_2nd_order(initial,grid_resolution,comm)
  initial => zero
  half_uu = periodic_2nd_order(initial,grid_resolution,comm)
  u_half = periodic_2nd_order(initial,grid_resolution,comm)
  do while(t<=t_final) !2nd-order Runge-Kutta: 
   dt = u%runge_kutta_2nd_step(nu ,grid_resolution)
    half_uu = u*u*half
    u_half = u + (u%xx()*nu - half_uu%x())*dt*half ! first substep
    half_uu = u_half*u_half*half
    u  = u + (u_half%xx()*nu - half_uu%x())*dt ! second substep
    t = t + dt
    if (comm%MyPID()==0) print *,'timestep=',tstep,'time=',t
    tstep=tstep+1
  end do
  if (comm%MyPID()==0) print *,'u at t=',t
#ifdef HAVE_MPI
  t_end= MPI_Wtime()
#else
  call cpu_time(t_end)
#endif
  if (comm%MyPID()==0) print *,'Elapsed CPU time=',t_end-t_start
  call u%output(comm)
  call half_uu%force_finalize
  call u_half%force_finalize
  call u%force_finalize
  call comm%force_finalize
#ifdef HAVE_MPI
  call MPI_FINALIZE(rc)
#endif
end program
