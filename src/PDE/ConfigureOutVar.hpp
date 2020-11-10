// *****************************************************************************
/*!
  \file      src/PDE/ConfigureOutVar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Assign functions to compute output variables from numerical
               solution
  \details   Assign functions to compute output variables from numerical
               solution.
*/
// *****************************************************************************
#ifndef ConfigureOutVar_h
#define ConfigureOutVar_h

#include <vector>

#include "FunctionPrototypes.hpp"

namespace inciter {

//! Assign function that computes output variables from numerical solution
tk::GetVarFn
assignGetVar( const std::string& name );

} // inciter::

#endif // ConfigureOutVar_h
