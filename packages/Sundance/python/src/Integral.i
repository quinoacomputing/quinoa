// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceIntegral.hpp"
#include "SundanceEssentialBC.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

%rename(Integral) makeIntegral;
%rename(EssentialBC) makeEssentialBC;


%inline %{
  /* */
  Sundance::Expr makeIntegral(const Sundance::CellFilter& domain,
                                  const Sundance::Expr& integrand,
                                  const Sundance::QuadratureFamily& quad)
  {
    return Sundance::Integral(domain, integrand, quad);
  }
  %}



%inline %{
  /* */
  Sundance::Expr makeEssentialBC(const Sundance::CellFilter& domain,
                                     const Sundance::Expr& integrand,
                                     const Sundance::QuadratureFamily& quad)
  {
    return Sundance::EssentialBC(domain, integrand, quad);
  }
  %}


