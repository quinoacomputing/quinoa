module vector_field_module
#include "ForTrilinos_config.h"
  use ForTrilinos_assertion_utility ,only : error_message,assert
  use FEpetra_Comm  ,only: Epetra_Comm
  use field_module  ,only : initial_field
  use scalar_field_module ,only: scalar_field
  use iso_c_binding       ,only : c_double
  use globals             ,only : nspace
  implicit none
  private
  public :: vector_field
  type  :: vector_field
    private
    type(scalar_field) ,allocatable ,dimension(:) :: u
    logical :: constructed=.false.
  contains
    procedure :: isConstructed
    procedure :: total
    procedure :: difference
    procedure :: negative
    procedure :: multiple
    procedure :: laplacian
    procedure :: dotGradient
    generic   :: operator(+) => total
    generic   :: operator(-) => difference, negative
    generic   :: operator(*) => multiple
    generic   :: operator(.laplacian.) => laplacian
    generic   :: operator(.dotGradient.) => dotGradient
    procedure :: runge_kutta_stable_step => rk_dt
    procedure :: output
    procedure :: force_finalize
  end type

  class(Epetra_Comm) , allocatable :: comm
  
 interface vector_field
    procedure new_vector_field ,default_vector_field
  end interface 

contains
  logical function isConstructed(this)
    class(vector_field) ,intent(in) :: this
    isConstructed = this%constructed
  end function

  function new_vector_field(initial_u,initial_v,initial_w,comm_in) result(this)
    type(vector_field)                       :: this
    procedure(initial_field)  ,pointer       :: initial_u,initial_v,initial_w
    class(Epetra_Comm)           ,intent(in) :: comm_in

    allocate(this%u(nspace))
    this%u(1) = scalar_field(initial_u,comm_in)
    this%u(2) = scalar_field(initial_v,comm_in)
    this%u(3) = scalar_field(initial_w,comm_in)

    if (.not.allocated(comm)) allocate(comm,source=comm_in)

    this%constructed=.true.

    ! Ensures (postcondition):
    call assert(this%isConstructed(),error_message('vector_field (constructor): construction failed.'))
  end function

  ! Default initialization function for vector field component
  function default_vector_field() result(this)
    use initializer ,only : zero
    type(vector_field)                  :: this
    procedure(initial_field)  ,pointer  :: default_initial

    allocate(this%u(nspace))
    default_initial=>zero
    this%u(1) = scalar_field(default_initial,comm)
    this%u(2) = scalar_field(default_initial,comm)
    this%u(3) = scalar_field(default_initial,comm)

    this%constructed=.true.

    ! Ensures (postcondition):
    call assert(this%isConstructed(),error_message('vector_field (constructor): construction failed.'))
  end function

  real(c_double) function rk_dt(this,order,nu)
    class(vector_field) ,intent(in) :: this 
    real(c_double) ,intent(in) :: nu
    integer ,intent(in) :: order
    
    ! Requires (preconditions):
    call assert(nu>0.,error_message('vector_field%rk_dt: invalid viscosity.'))
   
    rk_dt = this%u(1)%runge_kutta_stable_step(order,nu)

  end function

  type(vector_field) function total(lhs,rhs)
    class(vector_field) ,intent(in) :: lhs, rhs

    ! Requires (preconditions):
    call assert(lhs%isConstructed(),error_message('vector_field%total(): unconstructed right-hand side.'))
    call assert(rhs%isConstructed(),error_message('vector_field%total(): unconstructed left-hand side.'))

    total = vector_field()
    total%u(1)=lhs%u(1)+rhs%u(1)
    total%u(2)=lhs%u(2)+rhs%u(2)
    total%u(3)=lhs%u(3)+rhs%u(3)
    total%constructed = .true.
    
    ! Ensures (postconditions):
    call assert(total%isConstructed(),error_message('vector_field%total(): unconstructed right-hand side.'))
  end function

  type(vector_field) function negative(rhs)
    class(vector_field) ,intent(in) :: rhs

    ! Requires (preconditions):
    call assert(rhs%isConstructed(),error_message('vector_field%negative(): unconstructed right-hand side.'))

    negative = vector_field()
    negative%u(1)=-rhs%u(1)
    negative%u(2)=-rhs%u(2)
    negative%u(3)=-rhs%u(3)
    negative%constructed = .true.

    ! Ensures (postconditions):
    call assert(negative%isConstructed(),error_message('vector_field%negative(): negation failed.'))
  end function

  type(vector_field) function laplacian(rhs)
    class(vector_field) ,intent(in) :: rhs

    ! Requires (preconditions):
    call assert(rhs%isConstructed(),error_message('vector_field%laplacian(): unconstructed right-hand side.'))

    laplacian = vector_field()
    laplacian%u(1) = .laplacian.rhs%u(1)
    laplacian%u(2) = .laplacian.rhs%u(2)
    laplacian%u(3) = .laplacian.rhs%u(3)
    laplacian%constructed = .true.

    ! Ensures (postconditions):
    call assert(laplacian%isConstructed(),error_message('vector_field%laplacian(): laplacian operator failed.'))
  end function
 
  type(vector_field) function dotGradient(lhs,rhs)
    class(vector_field) ,intent(in) :: lhs, rhs

    ! Requires (preconditions):
    call assert(rhs%isConstructed(),error_message('vector_field%dotGradient(): unconstructed right-hand side.'))

    dotGradient = vector_field()
    dotGradient%u(1) = lhs%u(1)*rhs%u(1)%x() + lhs%u(2)*rhs%u(1)%y() + lhs%u(3)*rhs%u(1)%z()
    dotGradient%u(2) = lhs%u(1)*rhs%u(2)%x() + lhs%u(2)*rhs%u(2)%y() + lhs%u(3)*rhs%u(2)%z()
    dotGradient%u(3) = lhs%u(1)*rhs%u(3)%x() + lhs%u(2)*rhs%u(3)%y() + lhs%u(3)*rhs%u(3)%z()
    dotGradient%constructed = .true.

    ! Ensures (postconditions):
    call assert(dotGradient%isConstructed(),error_message('vector_field%dotGradient(): gradient operator failed.'))
  end function

  type(vector_field) function difference(lhs,rhs)
    class(vector_field) ,intent(in) :: lhs,rhs

    ! Requires (preconditions):
    call assert(lhs%isConstructed(),error_message('vector_field%difference(): unconstructed left-hand side.'))
    call assert(rhs%isConstructed(),error_message('vector_field%difference(): unconstructed right-hand side.'))

    difference = vector_field()
    difference%u(1) = lhs%u(1)-rhs%u(1)
    difference%u(2) = lhs%u(2)-rhs%u(2)
    difference%u(3) = lhs%u(3)-rhs%u(3)
    difference%constructed = .true.

    ! Ensures (postconditions):
    call assert(difference%isConstructed(),error_message('vector_field%difference(): difference operator failed.'))
  end function

  type(vector_field) function multiple(lhs,rhs)
    class(vector_field) ,intent(in) :: lhs
    real(c_double)     ,intent(in)  :: rhs

    ! Requires (preconditions):
    call assert(lhs%isConstructed(),error_message('vector_field%multiple(): unconstructed left-hand side.'))

    multiple = vector_field()
    multiple%u(1) = lhs%u(1)*rhs
    multiple%u(2) = lhs%u(2)*rhs
    multiple%u(3) = lhs%u(3)*rhs
    multiple%constructed = .true.

    ! Ensures (postconditions):
    call assert(multiple%isConstructed(),error_message('vector_field%multiple(): multiple operator failed.'))
  end function

  subroutine output(this)
    class(vector_field) ,intent(in) :: this

    ! Requires (preconditions):
    call assert(this%isConstructed(),error_message('vector_field%output(): unconstructed argument.'))

    call this%u(1)%output() 
    call this%u(2)%output() 
    call this%u(3)%output() 
  end subroutine

  subroutine force_finalize(this)
    class(vector_field), intent(inout) :: this
    call this%u(1)%force_finalize
    call this%u(2)%force_finalize
    call this%u(3)%force_finalize
  end subroutine
   
end module
