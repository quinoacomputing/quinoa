module hermetic_interface
  private
  public :: hermetic
  type, abstract :: hermetic
  contains
      procedure(free_memory), deferred :: cpp_delete
  end type
  abstract interface
     subroutine free_memory (this)
        import :: hermetic
        class(hermetic), intent(inout) :: this
     end subroutine
  end interface
end module
module ref_counter_implementation
  use hermetic_interface ,only : hermetic
  private; public :: ref_counter
  type ref_counter
    private
    integer, pointer :: count => null()
    class(hermetic), pointer :: obj => null()
  contains
    procedure, non_overridable :: grab
    procedure, non_overridable :: release
    procedure :: assign
#ifndef ForTrilinos_DISABLE_FINAL_SUBROUTINES
    final :: finalize_ref_counter
#endif 
    generic :: assignment(=) => assign
  end type
  interface ref_counter
    module procedure new_ref_counter; end interface
contains
  function new_ref_counter(object)
    class(hermetic), intent(in) :: object
    type(ref_counter), allocatable :: new_ref_counter
    allocate (new_ref_counter); allocate (new_ref_counter%count, source=0)
    allocate (new_ref_counter%obj, source=object)
    call new_ref_counter%grab
 end function
  subroutine grab(this)
    class(ref_counter), intent(inout) :: this
    if (associated(this%count)) then
      this%count = this%count + 1
    else; stop 'Error in grab: count not associated'
    end if
   end subroutine
  subroutine release(this)
    class (ref_counter), intent(inout) :: this
    if (associated(this%count)) then
      this%count = this%count - 1
      if (this%count == 0) then
        call this%obj%cpp_delete; deallocate (this%count, this%obj)
      else; this%count => null(); this%obj => null()
      end if
    else; stop 'Error in release: count not associated'
    end if
  end subroutine
  subroutine assign (lhs, rhs)
    class (ref_counter), intent(inout) :: lhs
    class (ref_counter), intent(in) :: rhs
    lhs%count => rhs%count; lhs%obj => rhs%obj
    call lhs%grab
 end subroutine
  recursive subroutine finalize_ref_counter (this)
    type(ref_counter), intent(inout) :: this
    if (associated(this%count)) call this%release
 end subroutine
end module
module universal_interface
  use hermetic_interface ,only: hermetic
  use ref_counter_implementation ,only: ref_counter
  implicit none
  private
  public :: universal
  type ,abstract ,extends(hermetic) :: universal
    type(ref_counter) :: counter
  contains
    procedure, non_overridable :: force_finalize
    procedure, non_overridable :: register_self
  end type
contains
  subroutine force_finalize (this)
    class(universal), intent(inout) :: this
    call this%counter%release
  end subroutine
  subroutine register_self (this)
    class(universal), intent(inout) :: this
    this%counter = ref_counter(this)
  end subroutine
end module
module faux_cpp_server
  use iso_c_binding ,only: c_int
  implicit none
  private
  public :: cpp_new_vector
  public :: cpp_delete_vector
  public :: num_leaks
  enum, bind(C) 
    enumerator :: never_allocated,is_allocated,is_deallocated
  end enum
  integer(c_int), save :: unique_id=0_c_int
  integer, parameter :: max_objects=1000
  integer(c_int), dimension(max_objects) :: cpp_vector=never_allocated
  interface cpp_new_vector ! public constructor
    module procedure initialize_vector
  end interface
contains
  function uninitialized_vector_id() result(id) 
    integer(c_int) :: id
    unique_id = unique_id + 1_c_int
    id = unique_id
  end function
  function initialize_vector() result(id) bind(C)
    integer(c_int) :: id
    id = uninitialized_vector_id()
    if (id>max_objects) stop 'Out of memory.'
    cpp_vector(id)=is_allocated
  end function
  subroutine cpp_delete_vector(id) bind(C)
    integer(c_int),value :: id
    cpp_vector(id) = is_deallocated
  end subroutine
  integer(c_int) function num_leaks()
    integer(c_int) :: num_never_allocated,num_deallocated
    num_never_allocated = count( cpp_vector == never_allocated )
    num_deallocated = count( cpp_vector == is_deallocated ) 
    num_leaks = max_objects - (num_never_allocated + num_deallocated)
  end function
end module
module vector_implementation
  use iso_c_binding, only : c_int
  use universal_interface, only: universal
  use faux_cpp_server
  implicit none

  private           ! Hide everything by default
  public :: vector  ! Expose type, constructor, type-bound procedures

  ! Shadow object
  type, extends(universal) :: vector 
    private
    integer(c_int) :: id ! C++ object identification tag
  contains
    procedure :: cpp_delete => call_cpp_delete_vector
  end type

  ! Constructors
  interface vector
     module procedure new_vector,default_vector
  end interface

contains

  type(vector) function default_vector(id)
    integer(c_int),intent(in) :: id
    default_vector%id = id
    call default_vector%register_self
  end function

  type(vector) function new_vector() 
    new_vector = vector(cpp_new_vector())
  end function

  subroutine call_cpp_delete_vector(this)
    class(vector),intent(inout) :: this
    call cpp_delete_vector(this%id)
  end subroutine
end module
program main
  use iso_fortran_env ,only : error_unit ,output_unit
  use faux_cpp_server ,only : num_leaks
  use vector_implementation, only: vector
  type(vector), allocatable, dimension(:) :: x   

  allocate(x(1))
  x(1) = vector() ! Object construction
  call x(1)%force_finalize() ! Final clean-up

  if (num_leaks()==0) then
    write(output_unit,*) 
    write(output_unit,fmt='(a)') "End Result: TEST PASSED" ;
  else
    write(error_unit,*) 
    write(error_unit,fmt='(a)') "End Result: TEST FAILED" ;
  end if
end program
