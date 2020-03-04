// *****************************************************************************
/*!
  \file      src/DiffEq/Dissipation/ConfigureDissipation.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the dissipation SDE
  \details   Register and compile configuration on the dissipation SDE.
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

#include "ConfigureDissipation.hpp"
#include "Dissipation.hpp"
#include "DissipationCoeffPolicy.hpp"

#include "CoupledEq.hpp"

namespace walker {

void
registerDissipation( DiffEqFactory& f, std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register dissipation SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for equation
  using DissipationPolicies =
    tk::cartesian_product< InitPolicies, DissipationCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< DissipationPolicies >(
    registerDiffEq< Dissipation >( f, ctr::DiffEqType::DISSIPATION, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoDissipation( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the dissipation SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  using eq = tag::dissipation;

  auto c = ++cnt[ ctr::DiffEqType::DISSIPATION ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DISSIPATION ), "" );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, eq, tag::depvar >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< eq >(c) ) );
  nfo.emplace_back( "number of components", std::to_string(
    g_inputdeck.get< tag::component, eq >()[c] ) );

  coupledInfo< eq, tag::velocity, tag::velocity_id >
             ( c, "velocity", nfo );

  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, eq, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, eq, tag::coeffpolicy >()[c] ) );

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, eq, tag::rng >()[c] ) );
  nfo.emplace_back( "coeff C3", std::to_string(
    g_inputdeck.get< tag::param, eq, tag::c3 >().at(c) ) );
  nfo.emplace_back( "coeff C4", std::to_string(
    g_inputdeck.get< tag::param, eq, tag::c4 >().at(c) ) );
  nfo.emplace_back( "coeff COM1", std::to_string(
    g_inputdeck.get< tag::param, eq, tag::com1 >().at(c) ) );
  nfo.emplace_back( "coeff COM2", std::to_string(
    g_inputdeck.get< tag::param, eq, tag::com2 >().at(c) ) );

  return nfo;
}

}  // walker::
