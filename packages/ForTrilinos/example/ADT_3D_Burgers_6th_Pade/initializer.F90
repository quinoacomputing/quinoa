module initializer
  use field_module ,only : initial_field
  use iso_c_binding, only: c_int, c_double
  use ForTrilinos_assertion_utility ,only : assert,error_message
      
  implicit none

contains

  real(c_double) pure function zero(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    zero = 0.
  end function

  real(c_double) pure function u_1D_initial(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    u_1D_initial = 10._c_double*sin(x)
  end function

  real(c_double) pure function u_3D_initial(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    u_3D_initial = 10._c_double*sin(x)
  end function

  real(c_double) pure function v_3D_initial(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    v_3D_initial = 10._c_double*sin(y)
  end function

  real(c_double) pure function w_3D_initial(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    w_3D_initial = 10._c_double*sin(z)
  end function

  real(c_double) pure function u_TaylorGreen2D(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    u_TaylorGreen2D = sin(x)*cos(y) 
  end function
  real(c_double) pure function v_TaylorGreen2D(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    v_TaylorGreen2D = -cos(x)*sin(y) 
  end function
  real(c_double) pure function vor3_TaylorGreen2D(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    vor3_TaylorGreen2D = 2.0*sin(x)*sin(y) 
  end function
  real(c_double) pure function Scalar3D(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    Scalar3D = sin(x)*cos(y)*cos(z) 
  end function
  real(c_double) pure function Scalar3D_dx(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    Scalar3D_dx = cos(x)*cos(y)*cos(z) 
  end function
  real(c_double) pure function Scalar3D_dy(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    Scalar3D_dy = -sin(x)*sin(y)*cos(z) 
  end function
  real(c_double) pure function Scalar3D_dz(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    Scalar3D_dz = -sin(x)*cos(y)*sin(z) 
  end function
  real(c_double) pure function Scalar3D_dx2(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    Scalar3D_dx2 = -sin(x)*cos(y)*cos(z) 
  end function
  real(c_double) pure function Scalar3D_dy2(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    Scalar3D_dy2 = -sin(x)*cos(y)*cos(z) 
  end function
  real(c_double) pure function Scalar3D_dz2(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    Scalar3D_dz2 = -sin(x)*cos(y)*cos(z) 
  end function
  real(c_double) pure function Scalar3D_laplacian(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    Scalar3D_laplacian = (-sin(x)*cos(y)*cos(z)) + (-sin(x)*cos(y)*cos(z)) + (-sin(x)*cos(y)*cos(z)) 
  end function
  real(c_double) pure function delta(x,y,z)
    real(c_double) ,intent(in) :: x,y,z
    real(c_double) :: r
    r = sqrt(x**2 + y**2 + z**2)
    if (r<epsilon(x)) then 
     delta = 1.
    else
     delta = 0.
    end if
  end function
end module 
