// -*- c++ -*-

%module PySundance

%feature("autodoc");

%exception 
{
  try
    {
      $action
    }
  catch (std::exception& e)
    {
      PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(e.what()));
      return NULL;
    }
}

%{
#include "Sundance.hpp"
#include "SundancePathUtils.hpp"
  %}


%inline %{


  void skipTimingOutput() {SundanceGlobal::skipTimingOutput()=true;}

  %}


%include "std_string.i"

namespace Sundance
{
  bool passFailTest(double err, double tol);
  std::string searchForFile(const std::string& name);
}


%inline%{
  class PyOut
  {
  public:
    PyOut() {}

    void write(const std::string& s) 
      {
        Sundance::Out::os() << s;
      }
  };
  %}

%include Mesh.i

%include Utils.i

%include Array.i

%include ParameterList.i

%include CellFilter.i

%include Playa.i

%include Quadrature.i

%include Basis.i

%include Spectral.i

%include Symbolics.i

%include CoordinateSystem.i

%include Integral.i

%include LinearProblem.i

%include NonlinearProblem.i

%include LinearEigenproblem.i

%include Functional.i

%include Viz.i

%include Discrete.i

%include AToC.i  

%include SpectralSolver.i  
