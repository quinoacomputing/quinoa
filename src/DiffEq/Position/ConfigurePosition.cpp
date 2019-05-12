// *****************************************************************************
/*!
  \file      src/DiffEq/Position/ConfigurePosition.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the position SDE
  \details   Register and compile configuration on the position SDE.
*/
// *****************************************************************************

#include <set>
#include <map>
#include <vector>
#include <string>
#include <utility>

#include <brigand/algorithms/for_each.hpp>

#include "Tags.hpp"
#include "CartesianProduct.hpp"
#include "DiffEqFactory.hpp"
#include "Walker/Options/DiffEq.hpp"
#include "Walker/Options/InitPolicy.hpp"

#include "ConfigurePosition.hpp"
#include "Position.hpp"
#include "PositionCoeffPolicy.hpp"

#include "CoupledEq.hpp"

namespace walker {

void
registerPosition( DiffEqFactory& f, std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register position SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for equation
  using PositionPolicies =
    tk::cartesian_product< InitPolicies, PositionCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< PositionPolicies >(
    registerDiffEq< Position >( f, ctr::DiffEqType::POSITION, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoPosition( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt )
// *****************************************************************************
//  Return information on the position SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  using eq = tag::position;

  auto c = ++cnt[ ctr::DiffEqType::POSITION ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::POSITION ), "" );

  nfo.emplace_back( "kind", "deterministic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, eq, tag::depvar >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< eq >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component, eq >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  coupledInfo< eq, tag::velocity, tag::velocity_id >
             ( c, "velocity", nfo );

  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, eq, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, eq, tag::coeffpolicy >()[c] ) );
  auto solve = g_inputdeck.get< tag::param, eq, tag::solve >()[c];
  nfo.emplace_back( "solve for", ctr::Depvar().name( solve ) );

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, eq, tag::rng >()[c] ) );

  return nfo;
}

}  // walker::
