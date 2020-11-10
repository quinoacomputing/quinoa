// *****************************************************************************
/*!
  \file      src/PDE/ConfigureOutVar.cpp
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

#include "ConfigureOutVar.hpp"
#include "ConfigureCompFlow.hpp"

tk::GetVarFn
inciter::assignGetVar( const std::string& name )
// *****************************************************************************
// Assign function that computes output variable from numerical solution
//! \param[in] name Name of variable whose OutVar::GetVarFn is to be assigned
//! \return Function assigned to output variable
// *****************************************************************************
{
  tk::GetVarFn f;

  assignCompFlowOutVar( name, f );

  // At this point all human-readable outvars must have a getvar fn assigned
  if (!name.empty())
    ErrChk( f, "OutVar::getvar() not assigned for output variable: " + name );

  return f;
}
