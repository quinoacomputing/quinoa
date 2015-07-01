// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceCoordinateSystem.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace Sundance
{
  

  /* */
  class CoordinateSystem
  {
  public:
    /* */
    Expr jacobian() const ;

    /** */
    bool supportsMeshDimension(int dim) const ;

  };
}

%rename(CartesianCoordinateSystem) makeCartesianCoordinateSystem;
%rename(MeridionalCylindricalCoordinateSystem) makeMeridionalCylindricalCoordinateSystem;
%rename(RadialSphericalCoordinateSystem) makeRadialSphericalCoordinateSystem;

%inline %{
  /* Create a maximal cell filter */
  Sundance::CoordinateSystem makeCartesianCoordinateSystem()
  {
    return new Sundance::CartesianCoordinateSystem();
  }
  %}
 
%inline %{
  /* Create a maximal cell filter */
  Sundance::CoordinateSystem makeMeridionalCylindricalCoordinateSystem()
  {
    return new Sundance::MeridionalCylindricalCoordinateSystem();
  }
  %}
 
%inline %{
  /* Create a maximal cell filter */
  Sundance::CoordinateSystem makeRadialSphericalCoordinateSystem()
  {
    return new Sundance::RadialSphericalCoordinateSystem();
  }
  %}
 
