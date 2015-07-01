module globals       ! Type kind parameter
  use iso_c_binding ,only :c_int,c_double
  implicit none
  private
  public :: nspace,grid_resolution,pi,x_length,tolerance,MaximumIter,dx,dy,dz,nx,ny,nz,threshold
  integer(c_int) ,parameter :: nspace = 3
  integer(c_int) ,parameter :: grid_resolution = 16
  real(c_double) ,parameter :: pi=acos(-1._c_double)
  real(c_double) ,parameter :: x_length=2.0*pi,y_length=x_length,z_length=x_length
  real(c_double) ,parameter :: tolerance = 1.0E-9
  real(c_double) ,parameter :: threshold = 1.0E-8
  integer(c_int) ,parameter :: MaximumIter = 100 ! Maximum number of iterations for solver Ax=b
  integer(c_int) ,parameter :: nx=grid_resolution,ny=nx, nz=nx
  real(c_double) ,parameter :: dx=x_length/real(nx,c_double)
  real(c_double) ,parameter :: dy=y_length/real(ny,c_double)  !if dirichlet in y dy=y_length/real(ny-1,c_double)
  real(c_double) ,parameter :: dz=z_length/real(nz,c_double)
end module
