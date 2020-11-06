// *****************************************************************************
/*!
  \file      src/Inciter/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Extract field output for inciter
  \details   Extract field output for inciter.
*/
// *****************************************************************************

#include "FieldOutput.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Tags.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

std::vector< std::string >
userFieldNames()
// *****************************************************************************
// Collect field output names based on user input
//! \return Output field names requested by user
// *****************************************************************************
{
  // Get list of output variables requested by user
  const auto& outvar = g_inputdeck.get< tag::cmd, tag::io, tag::outvar >();

  std::vector< std::string > nf;
  for (const auto& v : outvar) {
    std::stringstream s;
    s << v;
    nf.push_back( s.str() );
  }

  return nf;
}

std::vector< std::vector< tk::real > >
userFieldOutput( const tk::Fields& Un )
// *****************************************************************************
// Collect field output from solution based on user input
//! \return Output fields requested by user
// *****************************************************************************
{
  // Get list of output variables requested by user
  const auto& outvar = g_inputdeck.get< tag::cmd, tag::io, tag::outvar >();

  // Get offset map
  const auto& offset =
    g_inputdeck.get< tag::component >().offsetmap( g_inputdeck );

  std::vector< std::vector< tk::real > > nf;
  for (const auto& v : outvar)
    nf.push_back( Un.extract( v.field, tk::cref_find(offset,v.var) ) );

  return nf;
}

} // inciter::
