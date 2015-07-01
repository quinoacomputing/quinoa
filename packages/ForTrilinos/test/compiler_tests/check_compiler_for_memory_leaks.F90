module hermetic_interface
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
#include "ForTrilinos_config.h"
  use hermetic_interface ,only : hermetic
  type ref_counter
    integer, pointer :: count => null()
    class(hermetic), pointer :: obj => null()
  contains
    procedure, non_overridable :: grab
    procedure, non_overridable :: release
    procedure :: assign
#ifndef ForTrilinos_DISABLE_FINAL_SUBROUTINES
    final :: finalize_ref_counter
#endif /* ForTrilinos_DISABLE_FINAL_SUBROUTINES */
    generic :: assignment(=) => assign
  end type
contains
  function new_ref_counter(object)
    class(hermetic), intent(in) :: object
    type(ref_counter), allocatable :: new_ref_counter
    allocate (new_ref_counter); allocate (new_ref_counter%count, source=0)
    allocate (new_ref_counter%obj, source=object)
    call new_ref_counter%grab; end function
  subroutine grab(this)
    class(ref_counter), intent(inout) :: this
    if (associated(this%count)) then
      this%count = this%count + 1
    else; stop 'Error in grab: count not associated'
    end if; end subroutine
  subroutine release(this)
    class (ref_counter), intent(inout) :: this
    if (associated(this%count)) then
      this%count = this%count - 1
      if (this%count == 0) then
        call this%obj%cpp_delete
        deallocate (this%count); this%count => null()
        !deallocate (this%obj); this%obj => null() !This is the behaviour we want but  we need a workaround to avoid having a nested recursive function
        this%obj => null() ! workaround to avoid nested recursive function (small memory leak)
      else; this%count => null(); this%obj => null()
      end if
    else; stop 'Error in release: count not associated'
    end if; end subroutine
  subroutine assign (lhs, rhs)
    class (ref_counter), intent(inout) :: lhs
    class (ref_counter), intent(in) :: rhs
    lhs%count => rhs%count; lhs%obj => rhs%obj
    call lhs%grab; end subroutine
  recursive subroutine finalize_ref_counter (this)
    type(ref_counter), intent(inout) :: this
    if (associated(this%count)) call this%release; end subroutine
end module

module universal_interface
  use hermetic_interface ,only: hermetic
  use ref_counter_implementation ,only: ref_counter,new_ref_counter
  implicit none
  type ,abstract ,extends(hermetic) :: universal
    type(ref_counter) :: counter
  contains
    procedure, non_overridable :: force_finalize
    procedure, non_overridable :: register_self
    procedure                  :: component_finalization
  end type
contains
  recursive subroutine force_finalize (this)
    class(universal), intent(inout) :: this
    call this%component_finalization
    call this%counter%release
  end subroutine
  subroutine component_finalization(this)
    class(universal), intent(inout) :: this
  end subroutine
  subroutine register_self (this)
    class(universal), intent(inout) :: this
    this%counter = new_ref_counter(this)
  end subroutine
end module

module faux_cpp_server
  use iso_c_binding ,only: c_int,c_double
  implicit none
  integer(c_int), save :: unique_id=0_c_int
  enum, bind(C)
    enumerator :: never_allocated,is_allocated,is_deallocated
  end enum
  integer, parameter ::     max_objects=1000
  integer(c_int), dimension(max_objects) :: cpp_object=never_allocated
contains
  function uninitialized_object_id() result(id)
    integer(c_int) :: id
    unique_id = unique_id + 1_c_int
    id = unique_id
  end function
  function cpp_initialize_object() result(id) bind(C)
    integer(c_int) :: id
    id = uninitialized_object_id()
    if (id>max_objects) stop 'Out of memory.'
    cpp_object(id)=is_allocated
  end function
  function copy_object(original_id) result(id) bind(C)
    integer(c_int), intent(in) :: original_id
    integer(c_int) :: id
    id = uninitialized_object_id()
    if (id>max_objects) stop 'Out of memory.'
    cpp_object(id) =  cpp_object(original_id)
  end function
  subroutine cpp_delete_object(id) bind(C)
    integer(c_int),value :: id
    cpp_object(id) = is_deallocated
  end subroutine
  integer(c_int) function num_currently_allocated()
    num_currently_allocated = count(cpp_object == is_allocated)
  end function
end module

module foo_component_implementation
  use faux_cpp_server, only: cpp_initialize_object, cpp_delete_object
  use universal_interface ,only: universal
  implicit none
  type, extends(universal) :: foo_component
    integer :: foo_component_id ! C++ object identification tag
  contains
     procedure :: cpp_delete => call_cpp_delete_foo_component
  end type
contains
  type(foo_component) function new_foo_component()
    new_foo_component%foo_component_id = cpp_initialize_object()
    call new_foo_component%register_self
  end function
  subroutine call_cpp_delete_foo_component(this)
    class(foo_component),intent(inout) :: this
    call cpp_delete_object(this%foo_component_id)
  end subroutine
end module

module foo_parent_implementation
  use faux_cpp_server, only : cpp_initialize_object, cpp_delete_object
  use universal_interface ,only: universal
  use foo_component_implementation ,only: foo_component,new_foo_component
  implicit none
  type, extends(universal) :: foo_parent
    type(foo_component) :: component
    integer :: foo_parent_id ! C++ object identification tag
  contains
     procedure :: component_finalization => component_finalization_foo_parent_component
     procedure :: cpp_delete => call_cpp_delete_foo_parent
  end type
contains
  type(foo_parent) function new_foo_parent()
    new_foo_parent%foo_parent_id = cpp_initialize_object()
    new_foo_parent%component = new_foo_component()
    call new_foo_parent%register_self
  end function
  subroutine component_finalization_foo_parent_component(this)
    class(foo_parent) ,intent(inout) :: this
    call this%component%force_finalize()
  end subroutine
  subroutine call_cpp_delete_foo_parent(this)
    class(foo_parent),intent(inout) :: this
    call cpp_delete_object(this%foo_parent_id)
  end subroutine
end module

module foo_implementation
  use faux_cpp_server, only: cpp_initialize_object, cpp_delete_object
  use foo_parent_implementation ,only: foo_parent,new_foo_parent
  implicit none
  type, extends(foo_parent) :: foo
    integer :: foo_id ! C++ object identification tag
  contains
     procedure :: cpp_delete => call_cpp_delete_foo
  end type
contains
  type(foo) function new_foo()
    new_foo%foo_id = cpp_initialize_object()
    new_foo%foo_parent = new_foo_parent()
    call new_foo%register_self
  end function
  subroutine call_cpp_delete_foo(this)
    class(foo),intent(inout) :: this
    call this%foo_parent%cpp_delete()
    call cpp_delete_object(this%foo_id)
  end subroutine
end module

program main
  use iso_fortran_env, only : output_unit,error_unit
  use faux_cpp_server,    only : num_currently_allocated
  use foo_implementation, only : foo,new_foo
  type(foo) :: object
  object = new_foo()
  call object%force_finalize
  if (num_currently_allocated()==0) then
    write(output_unit,*)
    write(output_unit,fmt='(a)') "End Result: TEST PASSED" ;
  else
    write(error_unit,*)
    write(error_unit,fmt='(a)') "End Result: TEST FAILED" ;
  end if
end program
