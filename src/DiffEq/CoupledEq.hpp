// *****************************************************************************
/*!
  \file      src/DiffEq/CoupledEq.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functionality for querying information on coupled equations
  \details   Functionality for querying information on coupled equations.
*/
// *****************************************************************************
#ifndef CoupledEq_h
#define CoupledEq_h

#include <limits>

namespace walker {

extern ctr::InputDeck g_inputdeck;

template< typename eq, typename coupledeq, typename id >
void coupledInfo( std::size_t system,
                  std::string&& coupled_eq_name,
                  std::vector< std::pair< std::string, std::string > > nfo )
// *****************************************************************************
//  Query and generate info on coupled equations
//! \tparam eq Tag of the equation that is coupled
//! \tparam coupledeq Tag of the equation that is coupled to equation 'eq'
//! \tparam id Tag to access the coupled equation 'eq' (relative) ids, see
//!   tk::grm::couple in Control/Walker/InputDeck/Grammar.h
//! \param[in] system Relative equation system id of equation 'eq'
//! \param[in] coupled_eq_name Coupled equation name
//! \param[in,out] nfo Info vector to augment
// *****************************************************************************
{
  static_assert( !std::is_same< eq, coupledeq >::value,
                 "Eq and coupled eq must differ" );

  const auto& coupled_eq =
    g_inputdeck.template get< tag::param, eq, coupledeq >();
  const auto& coupled_eq_id = g_inputdeck.template get< tag::param, eq, id >();

  Assert( coupled_eq.size() == coupled_eq_id.size(),
          "Size mismatch in coupled eq configuration" );

  nfo.emplace_back( "coupled " + coupled_eq_name + " depvar",
                    std::string( 1, coupled_eq[ system ] ) );

  if (coupled_eq[ system ] != '-' )
    nfo.emplace_back( "coupled " + coupled_eq_name + " id",
                      std::to_string( coupled_eq_id[ system ] ) );
}

template< typename eq, typename coupledeq >
bool coupled( std::size_t system )
// *****************************************************************************
//  Query if equation 'eq' is coupled to equation 'coupledeq'
//! \tparam eq Tag of the equation that is coupled
//! \tparam coupledeq Tag of the equation that is coupled to equation 'eq'
//! \param[in] system Relative equation system id of equation 'eq'
//! \return True if equation 'eq' is coupled to equation 'coupledeq'
// *****************************************************************************
{
  static_assert( !std::is_same< eq, coupledeq >::value,
                 "Eq and coupled eq must differ" );

  const auto& coupled_eq =
    g_inputdeck.template get< tag::param, eq, coupledeq >();

  if (coupled_eq[ system ] != '-' ) return true; else return false;
}

template< typename eq, typename coupledeq >
char depvar( std::size_t system )
// *****************************************************************************
//  Query dependent variable of equation 'eq' coupled to equation 'coupledq'
//! \tparam eq Tag of the equation that is coupled
//! \tparam coupledeq Tag of the equation that is coupled to equation 'eq'
//! \param[in] system Relative equation system id of equation 'eq'
//! \return Character (dependent variable) of equation coupled to equation 'eq'
//! \note If equation 'eq' is not coupled to equation 'coupledeq', '-' is
//!   returned
// *****************************************************************************
{
  static_assert( !std::is_same< eq, coupledeq >::value,
                 "Eq and coupled eq must differ" );

  const auto& coupled_eq =
    g_inputdeck.template get< tag::param, eq, coupledeq >();

  return coupled_eq[ system ];
}

template< typename eq, typename coupledeq, typename id >
std::size_t system_id( std::size_t system )
// *****************************************************************************
//  Query relative id of coupled equation of potentially multiple eqs
//! \tparam eq Tag of the equation that is coupled
//! \tparam coupledeq Tag of the equation that is coupled to equation 'eq'
//! \tparam id Tag to access the coupled equation 'eq' (relative) ids, see
//!   tk::grm::couple in Control/Walker/InputDeck/Grammar.h
//! \param[in] system Relative equation system id of equation 'eq'
//! \return Relative id of coupled equation of potentially multiple eqs
//! \note If equation 'eq' is not coupled to equation 'coupledeq', we return a
//!   large number which will hopefully trigger some problems in client code if
//!   used. This is used to indicate "no coupling" so that client code can still
//!   call this function even for an equation that is not coupled, without
//!   chaning client code, compared to equations that are coupled. In other
//!   words, calling this function on a coupledeq equation that is not coupled
//!   is not an error.
// *****************************************************************************
{
  static_assert( !std::is_same< eq, coupledeq >::value,
                 "Eq and coupled eq must differ" );

  const auto& coupled_eq_id = g_inputdeck.template get< tag::param, eq, id >();

  if (coupled< eq, coupledeq >( system ))
    return coupled_eq_id[ system ];
  else
    return std::numeric_limits< std::size_t >::max();
}

template< typename eq, typename coupledeq, typename id >
std::size_t offset( std::size_t system )
// *****************************************************************************
//  Query system offset of coupled equation in tk::Data array of all systems
//! \tparam eq Tag of the equation that is coupled
//! \tparam coupledeq Tag of the equation that is coupled to equation 'eq'
//! \tparam id Tag to access the coupled equation 'eq' (relative) ids, see
//!   tk::grm::couple in Control/Walker/InputDeck/Grammar.h
//! \param[in] system Relative equation system id of equation 'eq'
//! \return System offset of coupled equation in tk::Data array of all systems
//! \note If equation 'eq' is not coupled to equation 'coupledeq', we return the
//!   id of the coupled equation. In other words, calling this function on a
//!   coupledeq equation that is not coupled is not an error.
// *****************************************************************************
{
  static_assert( !std::is_same< eq, coupledeq >::value,
                 "Eq and coupled eq must differ" );

  // Query relative id of coupled eq
  auto cid = system_id< eq, coupledeq, id >( system );

  // Return offset for coupled eq
  if (coupled< eq, coupledeq >( system ))
    return g_inputdeck.get< tag::component >().offset< coupledeq >( cid );
  else 
    return cid;
}

template< typename eq, typename coupledeq, typename id >
std::size_t ncomp( std::size_t system )
// *****************************************************************************
//  Query number of components of coupled equation
//! \tparam eq Tag of the equation that is coupled
//! \tparam coupledeq Tag of the equation that is coupled to equation 'eq'
//! \tparam id Tag to access the coupled equation 'eq' (relative) ids, see
//!   tk::grm::couple in Control/Walker/InputDeck/Grammar.h
//! \param[in] system Relative equation system id of equation 'eq'
//! \return Number of scalar components of coupled equation
//! \note If equation 'eq' is not coupled to equation 'coupledeq', we return the
//!   id of the coupled equation. In other words, calling this function on a
//!   coupledeq equation that is not coupled is not an error.
// *****************************************************************************
{
  static_assert( !std::is_same< eq, coupledeq >::value,
                 "Eq and coupled eq must differ" );

  // Query relative id of coupled eq
  auto cid = system_id< eq, coupledeq, id >( system );

  // Return number of scalar components of coupled eq
  if (coupled< eq, coupledeq >( system ))
    return g_inputdeck.get< tag::component >().get< coupledeq >().at( cid );
  else
    return cid;
}

} // walker::

#endif // CoupledEq_h
