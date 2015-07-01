module field_module
  use iso_c_binding, only: c_double
  implicit none
  private
  public :: initial_field
  public :: field
  type ,abstract :: field
  end type  
  abstract interface
    real(c_double) pure function initial_field(x,y,z)
      import :: c_double
      real(c_double) ,intent(in) :: x,y,z
    end function
  end interface
end module
