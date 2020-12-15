// *****************************************************************************
/*!
  \file      src/PDE/ConfigureOutVar.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Assign functions to compute output variables from the numerical
               solution
  \details   Assign functions to compute output variables from the numerical
               solution.
*/
// *****************************************************************************

#include "ConfigureOutVar.hpp"
#include "ConfigureCompFlow.hpp"
#include "ConfigureMultiMat.hpp"
#include "ConfigureTransport.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

//! Function object for counting the total number of eq systems configured
struct Neq {
  std::size_t& neq;
  explicit Neq( std::size_t& n ) : neq(n) {}
  template< typename U > void operator()( brigand::type_<U> ) {
    neq += inciter::g_inputdeck.get< tag::param, U, tag::depvar >().size();
  }
};

tk::GetVarFn
inciter::assignGetVars( const std::string& name )
// *****************************************************************************
// Assign all functions that compute output variables from numerical solutions
//! \param[in] name Name of variable whose OutVar::GetVarFn is to be assigned
//! \return Function assigned to output variable
//! \note This is used to configure human-readable output variables only.
// *****************************************************************************
{
  tk::GetVarFn f;

  // Query total number of eq sytems configured by user
  std::size_t neq = 0;
  brigand::for_each< ctr::parameters::Keys >( Neq(neq) );

  // Only attempt to configure getvars if we are called after the inputdeck has
  // been populated. This guard is here because OutVars, stored in the
  // inputdeck, may also be called by the runtime system on an empty inputdeck,
  // in which case we do nothing, but wait for when we are called with the
  // inputdeck populated.
  if (neq) {
    if (!g_inputdeck.get< tag::param, tag::transport, tag::depvar >().empty())
      assignTransportGetVars( name, f );
    if (!g_inputdeck.get< tag::param, tag::compflow, tag::depvar >().empty())
      assignCompFlowGetVars( name, f );
    if (!g_inputdeck.get< tag::param, tag::multimat, tag::depvar >().empty())
      assignMultiMatGetVars( name, f );
    // At this point all non-analytic human-readable outvars must have a getvar
    // function assigned
    if (!name.empty() && name.find("analytic") == std::string::npos)
      ErrChk( f, "OutVar::getvar() not assigned for output variable: " + name );
  }

  return f;
}
