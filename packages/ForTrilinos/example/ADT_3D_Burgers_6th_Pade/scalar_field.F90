! Scalar field implementation for a spacial discretization 
! and differentialtion with 6th order Pade
module scalar_field_module
#include "ForTrilinos_config.h"
  use ForTrilinos_assertion_utility ,only: error_message,assert,assert_identical
  use FEpetra_Comm     ,only: Epetra_Comm
  use FEpetra_Map      ,only: Epetra_Map
  use FEpetra_Vector   ,only: Epetra_Vector
  use FEpetra_CrsMatrix,only: Epetra_CrsMatrix
  use FAztecOO         ,only: AztecOO
  use ForTrilinos_error,only: error
  use field_module  ,only : initial_field
  use iso_c_binding ,only : c_double,c_int
  use globals       ,only : nspace, tolerance, MaximumIter ,nx,ny,nz,dx
  implicit none
  private
  public :: scalar_field
  type   :: scalar_field
    private
    type(Epetra_Vector) :: f
    logical :: constructed=.false.
  contains
    procedure :: isConstructed
    procedure :: total
    procedure :: difference
    procedure :: negative
    procedure :: product
    procedure :: multiple
    procedure :: laplacian
    procedure :: gradient
    generic   :: operator(+) => total
    generic   :: operator(-) => difference, negative
    generic   :: operator(*) => product, multiple
    generic   :: operator(.laplacian.) => laplacian
    generic   :: operator(.gradient.) => gradient
    procedure :: runge_kutta_stable_step => rk_dt
    procedure :: x  => df_dx     ! 1st derivative w.r.t. x
    procedure :: xx => d2f_dx2   ! 2nd derivative w.r.t. x
    procedure :: y  => df_dy     ! 1st derivative w.r.t. y
    procedure :: yy => d2f_dy2   ! 2nd derivative w.r.t. y
    procedure :: z  => df_dz     ! 1st derivative w.r.t. z
    procedure :: zz => d2f_dz2   ! 2nd derivative w.r.t. z
    procedure :: output
    procedure :: abs_maximum 
    procedure :: force_finalize
  end type

  ! Module variables 
  type(Epetra_Map) ,allocatable :: map  
  real(c_double) ,dimension(nx) :: x_node
  integer(c_int) ,save :: NumMy_xy_planes,my_first_xy_plane,my_last_xy_plane

  ! Module constants 
  integer(c_int) ,parameter :: IndexBase=1
  integer(c_int) ,parameter :: NumGlobalElements = nx*ny*nz
  integer(c_int) ,parameter :: NumGlobal_xy_planes = nz, Num_xy_points_per_plane = nx*ny
  ! Coefficients from equation (2.1) in S. Lele, J. Comp. Phys. 103 (1992)
  ! coef_1st=[beta, alpha, alpha, beta, c, b, a] 1st derivative coefficients
  ! coef_2nd=[beta, alpha, alpha, beta, c, b, a] 2nd derivative coefficients
  integer ,parameter :: max_stencil_width=7,rhs_stencil_width=4
  real(c_double) ,dimension(max_stencil_width) ,parameter :: coef_1st=[0., 1./3. , 1./3. , 0., 0. ,1./9. , 14./9. ]
  real(c_double) ,dimension(max_stencil_width) ,parameter :: coef_2nd=[0., 2./11., 2./11., 0., 0. ,3./11., 12./11.]

  interface scalar_field
    procedure new_scalar_field
  end interface 

contains
  logical function isConstructed(this)
    class(scalar_field) ,intent(in) :: this
    isConstructed = this%constructed
  end function

  function new_scalar_field(initial,comm) result(this)
    type(scalar_field)                          :: this
    procedure(initial_field) ,pointer           :: initial
    class(Epetra_Comm)              ,intent(in) :: comm
    real(c_double) ,dimension(:) ,allocatable   :: f_v
    integer(c_int) :: i,j,k,NumMyElements
 
    ! Requires (preconditions):
    call assert(mod(NumGlobal_xy_planes,comm%NumProc())==0 &
      ,error_message('scalar_field (constructor): number of processes not divisible by number of xy planes.'))

    ! Define local variables:
    NumMy_xy_planes = NumGlobal_xy_planes/comm%NumProc()
    NumMyElements = NumMy_xy_planes*Num_xy_points_per_plane

    ! Define module variables:
    forall(i=1:nx) x_node(i) =(i-1)*dx
    if (.not. allocated(map)) then
      map = Epetra_Map(NumGlobalElements,NumMyElements,IndexBase,comm)
    end if

    ! Define derived type components:
    ! Initialize and spread 3D field values along linear 1D array 
    allocate(f_v(NumMyElements))
    my_first_xy_plane = comm%MyPID()*NumMy_xy_planes + 1  
    my_last_xy_plane = (comm%MyPID()+1)*NumMy_xy_planes
    forall(i=1:nx,j=1:ny,k=my_first_xy_plane:my_last_xy_plane) &
      f_v( i + (j-1)*ny + (k-my_first_xy_plane)*ny*nz)=initial(x_node(i),x_node(j),x_node(k)) 

    this%f=Epetra_Vector(map,zero_initial=.true.)
    call this%f%ReplaceGlobalValues(f_v,map%MyGlobalElements())

    this%constructed=.true.

    ! Ensures (postcondition):
    call assert(this%isConstructed(),error_message('scalar_field (constructor): construction failed.'))
  end function

  real(c_double) function rk_dt(this,order,nu)
    class(scalar_field) ,intent(in) :: this 
    real(c_double) ,intent(in) :: nu
    integer ,intent(in) :: order
    real(c_double) ,parameter :: k_max=nx/2.0_c_double,safety_factor=0.9_c_double
    ! See eq. (3.1.10) of  Lele (J. Comp. Phys., 1992) 
    real(c_double) :: viscous_param = &
      ( (1+2.0*coef_2nd(2)*cos(k_max*dx)+2.0*coef_2nd(1)*cos(2.0*k_max*dx))/ &
        (  (2.0*coef_2nd(7)*(1-cos(k_max*dx))        ) +      &
           (coef_2nd(6)*(1-cos(2.0*k_max*dx))/2.0    ) +      &
           (2.0*coef_2nd(5)*(1-cos(3.0*k_max*dx))/9.0) )  )
    real(c_double) :: stable_viscous_param

    ! Requires (preconditions):
    call assert(nu>0.,error_message('scalar_field%rk_dt: invalid viscosity.'))

    ! See eq. (4.6.4) of  Lele (J. Comp. Phys., 1992) 
     select case(order)
       case(2)
         stable_viscous_param = 2.*viscous_param
       case(3)
         stable_viscous_param = 2.9*viscous_param
       case default
         call final_stop('Runge-Kutta order not supported.')
     end select
    rk_dt = (safety_factor*stable_viscous_param*(dx**2))/(nu*real(nspace,c_double))
   
    ! Ensures (postconditions):
    call assert(rk_dt < (stable_viscous_param*(dx**2))/(nu*real(nspace,c_double)), &
      error_message('scalar_field%rk_dt: time step too large for RK3 marching.'))
  end function

  type(scalar_field) function total(lhs,rhs)
    class(scalar_field) ,intent(in) :: lhs,rhs
    real(c_double), parameter :: one=1._c_double, zero=0._c_double

    ! Requires (preconditions):
    call assert(lhs%isConstructed(),error_message('scalar_field%total(): unconstructed right-hand side.'))
    call assert(rhs%isConstructed(),error_message('scalar_field%total(): unconstructed left-hand side.'))

    total%f=Epetra_Vector(map,zero_initial=.true.)
    call total%f%Update(one,lhs%f,one,rhs%f,zero)
    total%constructed = .true.
    
    ! Ensures (postconditions):
    call assert(total%isConstructed(),error_message('scalar_field%total(): unconstructed right-hand side.'))
  end function

  type(scalar_field) function negative(rhs)
    class(scalar_field) ,intent(in) :: rhs
    real(c_double), parameter :: minus_one=-1._c_double

    ! Requires (preconditions):
    call assert(rhs%isConstructed(),error_message('scalar_field%negative(): unconstructed right-hand side.'))

    negative%f=Epetra_Vector(map,zero_initial=.true.)
    call negative%f%Scale(minus_one,rhs%f)
    negative%constructed = .true.

    ! Ensures (postconditions):
    call assert(negative%isConstructed(),error_message('scalar_field%negative(): negation failed.'))
  end function

  type(scalar_field) function laplacian(rhs)
    class(scalar_field) ,intent(in) :: rhs

    ! Requires (preconditions):
    call assert(rhs%isConstructed(),error_message('scalar_field%laplacian(): unconstructed right-hand side.'))

    laplacian%f=Epetra_Vector(map,zero_initial=.true.)
    laplacian = rhs%xx() + rhs%yy() + rhs%zz()
    laplacian%constructed = .true.

    ! Ensures (postconditions):
    call assert(laplacian%isConstructed(),error_message('scalar_field%laplacian(): laplacian operator failed.'))
  end function
 
  type(scalar_field) function gradient(rhs)
    class(scalar_field) ,intent(in) :: rhs

    ! Requires (preconditions):
    call assert(rhs%isConstructed(),error_message('scalar_field%gradient(): unconstructed right-hand side.'))

    gradient = rhs%x() + rhs%y() + rhs%z()
    gradient%constructed = .true.

    ! Ensures (postconditions):
    call assert(gradient%isConstructed(),error_message('scalar_field%gradient(): gradient operator failed.'))
  end function

  type(scalar_field) function difference(lhs,rhs)
    class(scalar_field) ,intent(in) :: lhs,rhs
    real(c_double) ,parameter :: one=1._c_double, zero=0._c_double, minus_one=-1._c_double

    ! Requires (preconditions):
    call assert(lhs%isConstructed(),error_message('scalar_field%difference(): unconstructed left-hand side.'))
    call assert(rhs%isConstructed(),error_message('scalar_field%difference(): unconstructed right-hand side.'))

    difference%f=Epetra_Vector(map,zero_initial=.true.)
    call difference%f%Update(one,lhs%f,minus_one,rhs%f,zero)
    difference%constructed = .true.

    ! Ensures (postconditions):
    call assert(difference%isConstructed(),error_message('scalar_field%difference(): difference operator failed.'))
  end function

  type(scalar_field) function product(lhs,rhs)
    class(scalar_field) ,intent(in) :: lhs,rhs

    ! Requires (preconditions):
    call assert(lhs%isConstructed(),error_message('scalar_field%product(): unconstructed left-hand side.'))
    call assert(rhs%isConstructed(),error_message('scalar_field%product(): unconstructed right-hand side.'))

    product%f=Epetra_Vector(map,zero_initial=.true.)
    call product%f%Multiply(1._c_double,lhs%f,rhs%f,0._c_double)
    product%constructed=.true.

    ! Ensures (postconditions):
    call assert(product%isConstructed(),error_message('scalar_field%product(): product operator failed.'))
  end function

  type(scalar_field) function multiple(lhs,rhs)
    class(scalar_field) ,intent(in) :: lhs
    real(c_double) ,intent(in)  :: rhs

    ! Requires (preconditions):
    call assert(lhs%isConstructed(),error_message('scalar_field%multiple(): unconstructed left-hand side.'))

    multiple%f=Epetra_Vector(map,zero_initial=.true.)
    call multiple%f%Scale(rhs,lhs%f)
    multiple%constructed = .true.

    ! Ensures (postconditions):
    call assert(multiple%isConstructed(),error_message('scalar_field%multiple(): multiple operator failed.'))
  end function

  type(Epetra_Vector) function RHS_diff_x(this,coef,order) 
    class(scalar_field)         ,intent(in) :: this
    real(c_double),dimension(:) ,intent(in) :: coef
    integer(c_int)              ,intent(in) :: order
    real(c_double),dimension(:),allocatable :: Bf_v,f_v
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int)                          :: i, NumMyElements
    
    MyGlobalElements=map%MyGlobalElements()
    NumMyElements=map%NumMyElements()
    
    f_v=this%f%ExtractCopy() 
    associate (Bf => RHS_diff_x)
       Bf=Epetra_Vector(map)
       allocate(Bf_v(Bf%MyLength()))
       
       select case (order)
         case(1)     ! RHS for 1st derivative
            do i=1,NumMyElements
               if (mod(MyGlobalElements(i)-1,nx)==0) then
                  Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i+2)-f_v(i+nx-2))+(coef(7)/(2.0*dx))*(f_v(i+1)-f_v(i+nx-1))
               else if(mod(MyGlobalElements(i),nx)==0) then
                  Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i-nx+2)-f_v(i-2))+(coef(7)/(2.0*dx))*(f_v(i-nx+1)-f_v(i-1))
               else if (mod(MyGlobalElements(i)-2,nx)==0) then
                  Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i+2)-f_v(i+nx-2))+(coef(7)/(2.0*dx))*(f_v(i+1)-f_v(i-1))
               else if(mod(MyGlobalElements(i)+1,nx)==0) then
                  Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i-nx+2)-f_v(i-2))+(coef(7)/(2.0*dx))*(f_v(i+1)-f_v(i-1))
               else
                  Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i+2)-f_v(i-2))+(coef(7)/(2.0*dx))*(f_v(i+1)-f_v(i-1))
               end if
            enddo
         case(2)   ! RHS for 2nd derivative
            do i=1,NumMyElements
               if (mod(MyGlobalElements(i)-1,nx)==0) then
                  Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i+2)-2.0*f_v(i)+f_v(i+nx-2))+ &
                            (coef(7)/(dx*dx)    )*(f_v(i+1)-2.0*f_v(i)+f_v(i+nx-1))
               else if(mod(MyGlobalElements(i),nx)==0) then
                  Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i-nx+2)-2.0*f_v(i)+f_v(i-2))+ &
                            (coef(7)/(dx*dx)    )*(f_v(i-nx+1)-2.0*f_v(i)+f_v(i-1))
               else if (mod(MyGlobalElements(i)-2,nx)==0) then
                  Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i+2)-2.0*f_v(i)+f_v(i+nx-2))+ &
                            (coef(7)/(dx*dx)    )*(f_v(i+1)-2.0*f_v(i)+f_v(i-1))
               else if(mod(MyGlobalElements(i)+1,nx)==0) then
                  Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i-nx+2)-2.0*f_v(i)+f_v(i-2))+ &
                            (coef(7)/(dx*dx)    )*(f_v(i+1)-2.0*f_v(i)+f_v(i-1))
               else
                  Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i+2)-2.0*f_v(i)+f_v(i-2))+ &
                            (coef(7)/(dx*dx)    )*(f_v(i+1)-2.0*f_v(i)+f_v(i-1))
               end if
             enddo
          case default
             stop 'RHS_diff_x: invalid derivative order.'
          end select
       call Bf%ReplaceGlobalValues(Bf_v,MyGlobalElements)
    end associate
  end function

  type(Epetra_CrsMatrix) function Matrix_diff_x(coef) 
    use ForTrilinos_enum_wrappers ,only: FT_Epetra_DataAccess_E_Copy
    real(c_double) ,dimension(:),intent(in) :: coef
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int),dimension(:),allocatable :: NumNz
    integer(c_int) :: NumMyElements,i
    integer(c_int) :: indices(4)
    real(c_double) :: values(4),one=1._c_double
    integer(c_int),parameter :: diagonal=1_c_int
    type(error) :: err

    ! Get update list and number of local equations from given Map
    NumMyElements = map%NumMyElements()
    call assert_identical( [NumGlobalElements,map%NumGlobalElements()] )
    MyGlobalElements = map%MyGlobalElements()

    ! Create an integer vector NumNz that is used to build the Epetra Matrix
    ! NumNz(i) is the number of non-zero elements for the ith global equation
    ! on this processor
    allocate(NumNz(NumMyElements))

    NumNz = 5

    associate (A => Matrix_diff_x)
      ! Create a Epetra_Matrix
      A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

      ! Add rows one at a time
      ! Need some vectors to help
      values(1:4)=coef(1:4)
      do i=1,NumMyElements
        if (mod(MyGlobalElements(i)-1,nx)==0) then
          indices(1) = MyGlobalElements(i)+2
          indices(2) = MyGlobalElements(i)+1
          indices(3) = MyGlobalElements(i)+nx-1
          indices(4) = MyGlobalElements(i)+nx-2 
        else if(mod(MyGlobalElements(i),nx)==0) then
          indices(1) = MyGlobalElements(i)-nx+2
          indices(2) = MyGlobalElements(i)-nx+1
          indices(3) = MyGlobalElements(i)-1
          indices(4) = MyGlobalElements(i)-2
        else if (mod(MyGlobalElements(i)-2,nx)==0) then
          indices(1) = MyGlobalElements(i)+2
          indices(2) = MyGlobalElements(i)+1
          indices(3) = MyGlobalElements(i)-1
          indices(4) = MyGlobalElements(i)+nx-2 
        else if(mod(MyGlobalElements(i)+1,nx)==0) then
          indices(1) = MyGlobalElements(i)-nx+2
          indices(2) = MyGlobalElements(i)+1
          indices(3) = MyGlobalElements(i)-1
          indices(4) = MyGlobalElements(i)-2
        else
          indices(1) = MyGlobalElements(i)+2
          indices(2) = MyGlobalElements(i)+1
          indices(3) = MyGlobalElements(i)-1
          indices(4) = MyGlobalElements(i)-2
        end if
        call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
        call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
        !Put in the diaogonal entry
        call A%InsertGlobalValues(MyGlobalElements(i),[one],[MyGlobalElements(i)],err)
        call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
      end do

      !Finish up
      call A%FillComplete(.true.,err)
    end associate
  end function

  type(scalar_field) function df_dx(this) 
    class(scalar_field)    ,intent(in)        :: this
    type(Epetra_CrsMatrix) ,allocatable       :: Ax
    type(Epetra_Vector)                       :: b
    type(AztecOO)                             :: Solver
    type(error)                               :: err
    integer(c_int) ,dimension(:) ,allocatable :: MyGlobalElements

    ! Requires (preconditions):
    call assert(this%isConstructed(),error_message('scalar_field%x(): unconstructed argument.'))

    ! Get update list and number of local equations from given Map
    MyGlobalElements = map%MyGlobalElements()
    
    associate ( x => df_dx%f )
       !create matrix
       if (.not.allocated(Ax)) Ax=Matrix_diff_x(coef_1st)
       !create vector x and b
       x=Epetra_Vector(Ax%RowMap())
       ! This should be replaced by a better guess, e.g. the value at the most recent substep
       call x%Random() 
       b=RHS_diff_x(this,coef_1st,1)

       Solver=AztecOO(Ax,x,b) 
       call Solver%iterate(Ax,x,b,MaximumIter,tolerance,err)
       call assert( [err%error_code()==0_c_int] , [error_message('Solver%iterate: failed')] )
    end associate
    df_dx%constructed = .true.
    call assert(same_type_as(df_dx,this),error_message('scalar_field%df_dx: type mismatch.'))

    ! Ensures (postconditions):
    call assert(df_dx%isConstructed(),error_message('scalar_field%x(): differentiation failed.'))
  end function

  type(scalar_field) function d2f_dx2(this) 
    class(scalar_field)           ,intent(in) :: this
    type(Epetra_Vector)                       :: b
    type(Epetra_CrsMatrix) ,allocatable       :: Axx
    type(AztecOO)                             :: Solver
    type(error)                               :: err
    integer(c_int) ,dimension(:) ,allocatable :: MyGlobalElements

    ! Requires (preconditions):
    call assert(this%isConstructed(),error_message('scalar_field%x(): unconstructed argument.'))

    ! Get update list and number of local equations from given Map
    MyGlobalElements = map%MyGlobalElements()

    associate ( x => d2f_dx2%f )
       !create matrix
       if (.not.allocated(Axx)) Axx=Matrix_diff_x(coef_2nd)

       !create vector x and  
       x=Epetra_Vector(Axx%RowMap())
       call x%Random()
       b=RHS_diff_x(this,coef_2nd,2)

       Solver=AztecOO(Axx,x,b) 
       call Solver%iterate(Axx,x,b,MaximumIter,tolerance,err)
       call assert( [err%error_code()==0_c_int] , [error_message('Solver%iterate: failed')] )
    end associate
    
    d2f_dx2%constructed = .true.
    call assert(same_type_as(d2f_dx2,this),error_message('scalar_field%d2f_dx2: type mismatch.'))

    ! Ensures (postconditions):
    call assert(d2f_dx2%isConstructed(),error_message('scalar_field%x(): differentiation failed.'))
  end function

  type(Epetra_Vector) function RHS_diff_y(this,coef,order) 
    class(scalar_field)         ,intent(in) :: this
    real(c_double) ,dimension(:),intent(in) :: coef
    integer(c_int)              ,intent(in) :: order
    class(Epetra_Comm)         ,allocatable :: comm
    real(c_double),dimension(:),allocatable :: Bf_v,f_v
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int) :: i,j,NumMyElements
    
    MyGlobalElements=map%MyGlobalElements()
    NumMyElements=map%NumMyElements()
    call map%Comm(comm)

    allocate(f_v(this%f%MyLength()))
    f_v=this%f%ExtractCopy() 
  
    associate ( Bf => RHS_diff_y )
       Bf=Epetra_Vector(map)
       allocate(Bf_v(Bf%MyLength()))

       select case (order) 
          case(1)   ! RHS for 1st derivative
            do j=1,nx/comm%NumProc()   ! loops through xy planes
              do i=(nx*nx*(j-1)+1),(nx*nx*(j-1)+nx)
                 Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i+2*nx)-f_v(i+nx*nx-2*nx))+ &
                           (coef(7)/(2.0*dx))*(f_v(i+nx)  -f_v(i+nx*nx-nx))
              enddo 
              do i=nx*nx*j-nx+1,nx*nx*j
                 Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i-nx*nx+2*nx)-f_v(i-2*nx))+ &
                           (coef(7)/(2.0*dx))*(f_v(i-nx*nx+nx)  -f_v(i-nx))
              enddo 
              do i=(nx*nx*(j-1)+1+nx),(nx*nx*(j-1)+2*nx)
                 Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i+2*nx)-f_v(i+nx*nx-2*nx))+ &
                           (coef(7)/(2.0*dx))*(f_v(i+nx)  -f_v(i-nx))
              enddo
              do i=nx*nx*j-2*nx+1,nx*nx*j-nx
                 Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i-nx*nx+2*nx)-f_v(i-2*nx))+ &
                           (coef(7)/(2.0*dx))*(f_v(i+nx)        -f_v(i-nx))
              enddo
              do i=nx*nx*(j-1)+2*nx+1,nx*nx*j-2*nx
                 Bf_v(i) = (coef(6)/(4.0*dx))*(f_v(i+2*nx)-f_v(i-2*nx))+ &
                           (coef(7)/(2.0*dx))*(f_v(i+nx)  -f_v(i-nx))
              enddo
            enddo
          case(2)  ! RHS for 2nd derivative
            do j=1,nx/comm%NumProc()   ! loops through xy planes
              do i=(nx*nx*(j-1)+1),(nx*nx*(j-1)+nx)
                 Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i+2*nx)-2.0*f_v(i)+f_v(i+nx*nx-2*nx))+ &
                           (coef(7)/(dx*dx)    )*(f_v(i+nx)  -2.0*f_v(i)+f_v(i+nx*nx-nx))
              enddo 
              do i=nx*nx*j-nx+1,nx*nx*j
                 Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i-nx*nx+2*nx)-2.0*f_v(i)+f_v(i-2*nx))+ &
                           (coef(7)/(dx*dx)    )*(f_v(i-nx*nx+nx)  -2.0*f_v(i)+f_v(i-nx))
              enddo 
              do i=(nx*nx*(j-1)+1+nx),(nx*nx*(j-1)+2*nx)
                 Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i+2*nx)-2.0*f_v(i)+f_v(i+nx*nx-2*nx))+ &
                           (coef(7)/(dx*dx)    )*(f_v(i+nx)  -2.0*f_v(i)+f_v(i-nx))
              enddo
              do i=nx*nx*j-2*nx+1,nx*nx*j-nx
                 Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i-nx*nx+2*nx)-2.0*f_v(i)+f_v(i-2*nx))+ &
                           (coef(7)/(dx*dx)    )*(f_v(i+nx)        -2.0*f_v(i)+f_v(i-nx))
              enddo
              do i=nx*nx*(j-1)+2*nx+1,nx*nx*j-2*nx
                 Bf_v(i) = (coef(6)/(4.0*dx*dx))*(f_v(i+2*nx)-2.0*f_v(i)+f_v(i-2*nx))+ &
                           (coef(7)/(dx*dx)    )*(f_v(i+nx)  -2.0*f_v(i)+f_v(i-nx))
              enddo 
            enddo
          case default 
            stop 'RHS_diff_y: invalid derivative order.'
       end select
       call Bf%ReplaceGlobalValues(Bf_v,MyGlobalElements)  
    end associate
  end function

  type(Epetra_CrsMatrix) function Matrix_diff_y(coef)
    use ForTrilinos_enum_wrappers ,only: FT_Epetra_DataAccess_E_Copy
    real(c_double) ,dimension(:), intent(in) :: coef
    class(Epetra_Comm)  ,allocatable         :: comm
    type(error) :: err
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int),dimension(:),allocatable :: NumNz
    integer(c_int) :: MyGlobalElements_diagonal(1)
    integer(c_int) :: NumMyElements,i,j,indices(4), NumEntries
    real(c_double) :: values(4),one=1._c_double
    integer(c_int),parameter :: diagonal=1_c_int

    ! Get update list and number of local equations from given Map
    NumMyElements = map%NumMyElements()
    call assert_identical( [NumGlobalElements,map%NumGlobalElements()] )
    MyGlobalElements = map%MyGlobalElements()
    call map%Comm(comm)
    ! Create an integer vector NumNz that is used to build the Epetra Matrix
    ! NumNz(i) is the number of non-zero elements for the ith global equation
    ! on this processor
    allocate(NumNz(NumMyElements))

    NumNz = 5

    associate ( A => Matrix_diff_y )
       ! Create a Epetra_Matrix
       A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

       ! Add rows one at a time
       ! Need some vectors to help
       values(1:4)=coef(1:4)
       do j=1,nx/comm%NumProc()   ! loops through xy planes
          do i=(nx*nx*(j-1)+1),(nx*nx*(j-1)+nx)
            indices(1)=MyGlobalElements(i)+2*nx
            indices(2)=MyGlobalElements(i)+nx
            indices(3)=MyGlobalElements(i)+nx*nx-nx
            indices(4)=MyGlobalElements(i)+nx*nx-2*nx
            call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
            call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
          enddo 
          do i=nx*nx*j-nx+1,nx*nx*j
            indices(1)=MyGlobalElements(i)-nx*nx+2*nx
            indices(2)=MyGlobalElements(i)-nx*nx+nx
            indices(3)=MyGlobalElements(i)-nx
            indices(4)=MyGlobalElements(i)-2*nx
            call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
            call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
          enddo 
          do i=(nx*nx*(j-1)+1+nx),(nx*nx*(j-1)+2*nx)
            indices(1)=MyGlobalElements(i)+2*nx
            indices(2)=MyGlobalElements(i)+nx
            indices(3)=MyGlobalElements(i)-nx
            indices(4)=MyGlobalElements(i)+nx*nx-2*nx
            call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
            call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
          enddo
          do i=nx*nx*j-2*nx+1,nx*nx*j-nx
            indices(1)=MyGlobalElements(i)-nx*nx+2*nx
            indices(2)=MyGlobalElements(i)+nx
            indices(3)=MyGlobalElements(i)-nx
            indices(4)=MyGlobalElements(i)-2*nx
            call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
            call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
          enddo
          do i=nx*nx*(j-1)+2*nx+1,nx*nx*j-2*nx
            indices(1) = MyGlobalElements(i)+2*nx
            indices(2) = MyGlobalElements(i)+nx
            indices(3) = MyGlobalElements(i)-nx
            indices(4) = MyGlobalElements(i)-2*nx
            call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
            call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
          enddo
        enddo
        do i=1,NumMyElements
          !Put in the diaogonal entry
          call A%InsertGlobalValues(MyGlobalElements(i),[one],[MyGlobalElements(i)],err)
          call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
        enddo

        !Finish up
        call A%FillComplete(.true.,err)
        call assert( [err%error_code()==0_c_int] , [error_message('A%FillComplete: failed')] )
     end associate
  end function

  type(scalar_field) function df_dy(this) 
    class(scalar_field)         ,intent(in) :: this
    type(Epetra_Vector)                     :: b
    type(Epetra_CrsMatrix),allocatable      :: Ay
    type(AztecOO)                           :: Solver
    type(error)                             :: err
    integer(c_int),dimension(:),allocatable :: MyGlobalElements

    ! Requires (preconditions):
    call assert(this%isConstructed(),error_message('scalar_field%y(): unconstructed argument.'))

    ! Get update list and number of local equations from given Map
    MyGlobalElements = map%MyGlobalElements()

    ! create matrix
    if (.not. allocated(Ay)) Ay=Matrix_diff_y(coef_1st)
 
    associate ( x => df_dy%f ) 
       !create vector x and b
       x=Epetra_Vector(Ay%RowMap())
       call x%Random()
       b=RHS_diff_y(this,coef_1st,1)

       Solver=AztecOO(Ay,x,b) 
       call Solver%iterate(Ay,x,b,MaximumIter,tolerance,err)
       call assert( [err%error_code()==0_c_int] , [error_message('Solver%iterate: failed')] )
    end associate
    df_dy%constructed = .true.

    ! Ensures (postconditions):
    call assert(df_dy%isConstructed(),error_message('scalar_field%y(): differentiation failed.'))
  end function
  
  type(scalar_field) function d2f_dy2(this) 
    class(scalar_field)         ,intent(in) :: this
    type(Epetra_Vector)                     :: b
    type(Epetra_CrsMatrix),allocatable      :: Ayy
    type(AztecOO)                           :: Solver
    type(error)                             :: err
    integer(c_int),dimension(:),allocatable :: MyGlobalElements

    ! Requires (preconditions):
    call assert(this%isConstructed(),error_message('scalar_field%yy(): unconstructed argument.'))

    ! Get update list and number of local equations from given Map
    MyGlobalElements = map%MyGlobalElements()

    ! create matrix
    if (.not. allocated(Ayy)) Ayy=Matrix_diff_y(coef_2nd)

    associate ( x => d2f_dy2%f )
       !create vector x and b  
       x=Epetra_Vector(Ayy%RowMap())
       call x%Random()
       b=RHS_diff_y(this,coef_2nd,2)

       Solver=AztecOO(Ayy,x,b) 
       call Solver%iterate(Ayy,x,b,MaximumIter,tolerance,err)
       call assert( [err%error_code()==0_c_int] , [error_message('Solver%iterate: failed')] )
    end associate
    d2f_dy2%constructed = .true.

    ! Ensures (postconditions):
    call assert(d2f_dy2%isConstructed(),error_message('scalar_field%yy(): differentiation failed.'))
  end function
  
  type(Epetra_Vector) function RHS_diff_z(this,coef,order) 
     use ForTrilinos_enum_wrappers,only: FT_Epetra_DataAccess_E_Copy
    class(scalar_field)          ,intent(in) :: this
    real(c_double) ,dimension(:) ,intent(in) :: coef
    integer(c_int)               ,intent(in) :: order
    type(Epetra_CrsMatrix)                   :: B
    type(error)                              :: err
    integer(c_int),dimension(:),allocatable  :: MyGlobalElements,NumNz
    real(c_double) :: diagonal_value, values(6)
    integer(c_int) :: indices(6),NumMyElements,i
   
    ! Creating target map and vector
    call assert_identical( [NumGlobalElements,map%NumGlobalElements()] )
    MyGlobalElements=map%MyGlobalElements()
    NumMyElements=map%NumMyElements()


    ! Create an integer vector NumNz that is used to build the Epetra Matrix
    ! NumNz(i) is the number of non-zero elements for the ith global equation
    ! on this processor
    allocate(NumNz(NumMyElements))

    NumNz = 7

    ! Create a Epetra_Matrix
    B = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

    ! Add rows one at a time
    ! Need some vectors to help
    select case (order) 
      case(1) ! RHS for 1st derivative
         values=(/coef(5)/(6.0*dx), coef(6)/(4.0*dx), coef(7)/(2.0*dx), &
               -coef(7)/(2.0*dx),-coef(6)/(4.0*dx),-coef(5)/(6.0*dx) /)
         diagonal_value = 0.0
      case(2) ! RHS for 2nd derivative
         values=(/coef(5)/(9.0*dx*dx), coef(6)/(4.0*dx*dx), coef(7)/(dx*dx), &
             coef(7)/(dx*dx)    , coef(6)/(4.0*dx*dx), coef(5)/(9.0*dx*dx) /)
         diagonal_value = -2.0*((coef(5)/(9.0*dx*dx))+(coef(6)/(4.0*dx*dx))+(coef(7)/(dx*dx)))
      case default 
         call final_stop('scalar_field%RHS_diff_z: This is not an implemented derivative')
    end select
    do i=1,NumMyElements
     if (MyGlobalElements(i).GE.1.and.MyGlobalElements(i).LE.nx*ny) then
        indices(1)=MyGlobalElements(i)+3*nx*ny
        indices(2)=MyGlobalElements(i)+2*nx*ny
        indices(3)=MyGlobalElements(i)+nx*ny
        indices(4)=MyGlobalElements(i)+nx*ny*nz-nx*ny
        indices(5)=MyGlobalElements(i)+nx*ny*nz-2*nx*ny
        indices(6)=MyGlobalElements(i)+nx*ny*nz-3*nx*ny
      else if (MyGlobalElements(i).GE.(nx*ny*(nz-1)+1).and.MyGlobalElements(i).LE.nx*ny*nz) then
        indices(1)=MyGlobalElements(i)-nx*ny*nz+3*nx*ny
        indices(2)=MyGlobalElements(i)-nx*ny*nz+2*nx*ny
        indices(3)=MyGlobalElements(i)-nx*ny*nz+nx*ny
        indices(4)=MyGlobalElements(i)-nx*ny
        indices(5)=MyGlobalElements(i)-2*nx*ny
        indices(6)=MyGlobalElements(i)-3*nx*ny
      else if (MyGlobalElements(i).GE.(nx*ny+1).and.MyGlobalElements(i).LE.2*nx*ny) then
        indices(1)=MyGlobalElements(i)+3*nx*ny
        indices(2)=MyGlobalElements(i)+2*nx*ny
        indices(3)=MyGlobalElements(i)+nx*ny
        indices(4)=MyGlobalElements(i)-nx*ny
        indices(5)=MyGlobalElements(i)+nx*ny*nz-2*nx*ny
        indices(6)=MyGlobalElements(i)+nx*ny*nz-3*nx*ny
      else if (MyGlobalElements(i).GE.(nx*ny*(nz-2)+1).and.MyGlobalElements(i).LE.(nx*ny*(nz-1))) then
        indices(1)=MyGlobalElements(i)-nx*ny*nz+3*nx*ny
        indices(2)=MyGlobalElements(i)-nx*ny*nz+2*nx*ny
        indices(3)=MyGlobalElements(i)+nx*ny
        indices(4)=MyGlobalElements(i)-nx*ny
        indices(5)=MyGlobalElements(i)-2*nx*ny
        indices(6)=MyGlobalElements(i)-3*nx*ny
      else if (MyGlobalElements(i).GE.(2*nx*ny+1).and.MyGlobalElements(i).LE.3*nx*ny) then
        indices(1)=MyGlobalElements(i)+3*nx*ny
        indices(2)=MyGlobalElements(i)+2*nx*ny
        indices(3)=MyGlobalElements(i)+nx*ny
        indices(4)=MyGlobalElements(i)-nx*ny
        indices(5)=MyGlobalElements(i)-2*nx*ny
        indices(6)=MyGlobalElements(i)+nx*ny*nz-3*nx*ny
      else if (MyGlobalElements(i).GE.(nx*ny*(nz-3)+1).and.MyGlobalElements(i).LE.(nx*ny*(nz-2))) then
        indices(1)=MyGlobalElements(i)-nx*ny*nz+3*nx*ny
        indices(2)=MyGlobalElements(i)+2*nx*ny
        indices(3)=MyGlobalElements(i)+nx*ny
        indices(4)=MyGlobalElements(i)-nx*ny
        indices(5)=MyGlobalElements(i)-2*nx*ny
        indices(6)=MyGlobalElements(i)-3*nx*ny
      else if (MyGlobalElements(i).GE.(3*nx*ny+1).and.MyGlobalElements(i).LE.(nx*ny*(nz-3))) then
        indices(1) = MyGlobalElements(i)+3*nx*ny
        indices(2) = MyGlobalElements(i)+2*nx*ny
        indices(3) = MyGlobalElements(i)+nx*ny
        indices(4) = MyGlobalElements(i)-nx*ny
        indices(5) = MyGlobalElements(i)-2*nx*ny
        indices(6) = MyGlobalElements(i)-3*nx*ny
      endif
      call B%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
      call assert( [err%error_code()==0_c_int] , [error_message('B%InsertGlobalValues: failed')] )
      !Put in the diaogonal entry
      call B%InsertGlobalValues(MyGlobalElements(i),[diagonal_value],[MyGlobalElements(i)],err)
      call assert( [err%error_code()==0_c_int] , [error_message('B%InsertGlobalValues: failed')] )
    enddo

    !Finish up
    call B%FillComplete(.true.,err)
    call assert( [err%error_code()==0_c_int] , [error_message('B%FillComplete: failed')] )
 
    associate ( Bf => RHS_diff_z )
       Bf=Epetra_Vector(map,zero_initial=.true.)
       call B%Multiply_Vector(.false.,this%f,Bf)
    end associate
  end function RHS_diff_z

  type(Epetra_CrsMatrix) function Matrix_diff_z(coef) 
    use ForTrilinos_enum_wrappers ,only: FT_Epetra_DataAccess_E_Copy
    real(c_double) ,dimension(:) ,intent(in) :: coef
    type(error)    :: err
    integer(c_int) :: NumMyElements,i,indices(4), NumEntries
    real(c_double) :: values(4),one=1._c_double
    integer(c_int) :: MyGlobalElements_diagonal(1)
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int),dimension(:),allocatable :: NumNz
    integer(c_int),parameter :: diagonal=1_c_int

    ! Get update list and number of local equations from given Map
    call assert_identical( [NumGlobalElements,map%NumGlobalElements()] )
    NumMyElements = map%NumMyElements()
    MyGlobalElements = map%MyGlobalElements()

    ! Create an integer vector NumNz that is used to build the Epetra Matrix
    ! NumNz(i) is the number of non-zero elements for the ith global equation
    ! on this processor
    allocate(NumNz(NumMyElements))

    NumNz = 5

    associate ( A => Matrix_diff_z) 
      ! Create a Epetra_Matrix
      A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

      ! Add rows one at a time
      ! Need some vectors to help
      values(1:4)=coef(1:4)
      do i=1,NumMyElements
        if (MyGlobalElements(i).GE.1.and.MyGlobalElements(i).LE.nx*nx) then
          indices(1)=MyGlobalElements(i)+2*nx*nx
          indices(2)=MyGlobalElements(i)+nx*nx
          indices(3)=MyGlobalElements(i)+nx*nx*nx-nx*nx
          indices(4)=MyGlobalElements(i)+nx*nx*nx-2*nx*nx
        else if (MyGlobalElements(i).GE.(nx*nx*(nx-1)+1).and.MyGlobalElements(i).LE.nx*nx*nx) then
          indices(1)=MyGlobalElements(i)-nx*nx*nx+2*nx*nx
          indices(2)=MyGlobalElements(i)-nx*nx*nx+nx*nx
          indices(3)=MyGlobalElements(i)-nx*nx
          indices(4)=MyGlobalElements(i)-2*nx*nx
        else if (MyGlobalElements(i).GE.(nx*nx+1).and.MyGlobalElements(i).LE.2*nx*nx) then
          indices(1)=MyGlobalElements(i)+2*nx*nx
          indices(2)=MyGlobalElements(i)+nx*nx
          indices(3)=MyGlobalElements(i)-nx*nx
          indices(4)=MyGlobalElements(i)+nx*nx*nx-2*nx*nx
        else if (MyGlobalElements(i).GE.(nx*nx*(nx-2)+1).and.MyGlobalElements(i).LE.(nx*nx*(nx-1))) then
          indices(1)=MyGlobalElements(i)-nx*nx*nx+2*nx*nx
          indices(2)=MyGlobalElements(i)+nx*nx
          indices(3)=MyGlobalElements(i)-nx*nx
          indices(4)=MyGlobalElements(i)-2*nx*nx
        else if (MyGlobalElements(i).GE.(2*nx*nx+1).and.MyGlobalElements(i).LE.(nx*nx*(nx-2))) then
          indices(1) = MyGlobalElements(i)+2*nx*nx
          indices(2) = MyGlobalElements(i)+nx*nx
          indices(3) = MyGlobalElements(i)-nx*nx
          indices(4) = MyGlobalElements(i)-2*nx*nx
        endif
        call A%InsertGlobalValues(MyGlobalElements(i),values,indices,err)
        call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
        !Put in the diaogonal entry
        call A%InsertGlobalValues(MyGlobalElements(i),[one],[MyGlobalElements(i)],err)
        call assert( [err%error_code()==0_c_int] , [error_message('A%InsertGlobalValues: failed')] )
      enddo

      !Finish up
      call A%FillComplete(.true.,err)
      call assert( [err%error_code()==0_c_int] , [error_message('A%FillComplete: failed')] )
    end associate
  end function

  type(scalar_field) function df_dz(this) 
    class(scalar_field)         ,intent(in) :: this
    type(Epetra_Vector)                     :: b
    type(Epetra_CrsMatrix),allocatable      :: Az
    type(AztecOO)                           :: Solver
    type(error)                             :: err
    integer(c_int),dimension(:),allocatable :: MyGlobalElements

    ! Requires (preconditions):
    call assert(this%isConstructed(),error_message('scalar_field%z(): unconstructed argument.'))

    ! Get update list and number of local equations from given Map
    MyGlobalElements = map%MyGlobalElements()

    ! Create a Epetra_Matrix
    if (.not.allocated(Az)) Az=Matrix_diff_z(coef_1st)
   
    associate ( x => df_dz%f ) 
      !create vector x and b
      x=Epetra_Vector(Az%RowMap())
      call x%Random()
      b=RHS_diff_z(this,coef_1st,1)
    
      Solver=AztecOO(Az,x,b) 
      call Solver%iterate(Az,x,b,MaximumIter,tolerance,err)
      call assert( [err%error_code()==0_c_int] , [error_message('Solver%iterate: failed')] )
    end associate
    df_dz%constructed = .true.

    ! Ensures (postconditions):
    call assert(df_dz%isConstructed(),error_message('scalar_field%z(): differentiation failed.'))
  end function
  
  type(scalar_field) function d2f_dz2(this) 
    class(scalar_field)         ,intent(in) :: this
    type(Epetra_Vector)                     :: b
    type(Epetra_CrsMatrix),allocatable      :: Azz
    type(AztecOO)                           :: Solver
    type(error)                             :: err
    integer(c_int),dimension(:),allocatable :: MyGlobalElements

    ! Requires (preconditions):
    call assert(this%isConstructed(),error_message('scalar_field%zz(): unconstructed argument.'))

    ! Get update list and number of local equations from given Map
    MyGlobalElements = map%MyGlobalElements()

    ! create matrix
    if (.not.allocated(Azz)) Azz=Matrix_diff_z(coef_2nd)

    associate (x => d2f_dz2%f )
      !create vector x and  
      x=Epetra_Vector(Azz%RowMap())
      call x%Random()
      b=RHS_diff_z(this,coef_2nd,2)

      Solver=AztecOO(Azz,x,b) 
      call Solver%iterate(Azz,x,b,MaximumIter,tolerance,err)
      call assert( [err%error_code()==0_c_int] , [error_message('Solver%iterate: failed')] )
    end associate
    d2f_dz2%constructed = .true.

    ! Ensures (postconditions):
    call assert(d2f_dz2%isConstructed(),error_message('scalar_field%zz(): differentiation failed.'))
  end function

  subroutine output(this)
    class(scalar_field) ,intent(in) :: this
    class(Epetra_Comm)  ,allocatable :: comm
    integer(c_int) :: NumMyElements
    integer(c_int) :: i,j,k,element
    integer(c_int), dimension(:), allocatable :: MyGlobalElements
    real(c_double), dimension(:), allocatable :: f_v
    real(c_double), dimension(:), allocatable :: f
    integer :: file_unit=20

    ! Requires (preconditions):
    call assert(this%isConstructed(),error_message('scalar_field%output(): unconstructed argument.'))

    call assert_identical( [NumGlobalElements,map%NumGlobalElements()] )
    NumMyElements=map%NumMyElements()
    MyGlobalElements=map%MyGlobalElements()
    call map%Comm(comm)

    allocate(f_v(this%f%MyLength()))
    f_v=this%f%ExtractCopy()
    allocate(f(NumGlobalElements))

    call comm%GatherAll(f_v,f)
    element=1
    do k=1,size(x_node)
     do j=1,size(x_node)
      do i=1,size(x_node)
        if (comm%MyPID()==0) write(file_unit,'(5(E20.12,1x))') x_node(i),x_node(j),x_node(k),f(element)
        element=element+1
      enddo
     enddo
    enddo
    file_unit=file_unit+1
    call assert([NumGlobalElements==element-1],[error_message('Error during output incorrect number of elements in Epetra_Vector')])
    call check_solution(f,NumGlobalElements,file_unit)
  end subroutine

  subroutine check_solution(f,NumGlobalElements,file_unit)
    integer(c_int) :: i,NumGlobalElements
    integer, intent(in) :: file_unit
    real(c_double), dimension(:), allocatable :: f,exact
    real(c_double),parameter :: tolerance=10.0E-1
    integer(c_int),parameter :: exact_grid=16
    logical :: success=.true.
    allocate(exact(exact_grid))
    exact=(/0.0, 0.686, 1.37, 2.053, 2.73, 3.387, 3.918, 3.57, & !exact solution at t=0.4626377
            0.0,-3.57,-3.918,-3.387,-2.73,-2.053,-1.37,-0.686/) !exact solution at t=0.4626377
    if (NumGlobalElements==exact_grid .and. file_unit==20) then
     if  (sqrt(sum((f(1:exact_grid)-exact(:))**2))>tolerance) success=.false.
    endif
    if (success) then
      print *
      print *, "End Result: TEST PASSED"
    else
      print *
      print *, "End Result: TEST FAILED"
    end if
  end subroutine 

  function abs_maximum(this) result(f)
   class(scalar_field), intent(in) :: this
   type(Epetra_Vector) :: f_abs
   real(c_double),allocatable,dimension(:) :: f
   f_abs=Epetra_Vector(map,zero_initial=.true.)
   call f_abs%Abs(this%f)
   f=f_abs%MaxValue()
   call assert( size(f)==1_c_int , error_message('MaxValue failed: return more than one value') )
  end function
  
  subroutine force_finalize(this)
    class(scalar_field), intent(inout) :: this
    call this%f%force_finalize
  end subroutine
   
  subroutine final_stop(message)
    use iso_fortran_env ,only : error_unit
    character(len=*) ,intent(in) :: message
    integer :: ierr
#ifdef HAVE_MPI
    call MPI_FINALIZE(ierr)
#endif
    write(error_unit,*)  message
    stop 'scalar_field (final_stop): message printed to error_unit.'
  end subroutine
end module
