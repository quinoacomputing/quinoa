// *****************************************************************************
/*!
  \file      src/PDE/ConfigureOutVar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Assign functions to compute output variables from the numerical
               solution
  \details   Assign functions to compute output variables from the numerical
               solution.
*/
// *****************************************************************************
#ifndef ConfigureOutVar_h
#define ConfigureOutVar_h

#include <vector>

#include "FunctionPrototypes.hpp"

namespace inciter {

//! \brief Assign all functions that compute output variables from the
//!   numerical solution
tk::GetVarFn
assignGetVars( const std::string& name );

//! Assign a function to compute an output variable from the numerical solution
//! \tparam Keyword Keyword used to match to variable name whose fn to assign
//! \param[in] name Name of variable whose OutVar::GetVarFn is to be assigned
//! \param[in] src Function to assign if there is a match
//! \param[in,out] dst Function to assign to if there is a match
//! \note This is used to configure human-readable output variables only.
template< class Keyword >
void
assign( const std::string& name, const tk::GetVarFn& src, tk::GetVarFn& dst ) {
  if (name == Keyword::string()) dst = src;
}

} // inciter::

#endif // ConfigureOutVar_h
