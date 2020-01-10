// *****************************************************************************
/*!
  \file      src/DiffEq/WrightFisher/ConfigureWrightFisher.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the Wright-Fisher SDE
  \details   Register and compile configuration on the Wright-Fisher SDE.
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

#include "ConfigureWrightFisher.hpp"
#include "WrightFisher.hpp"
#include "WrightFisherCoeffPolicy.hpp"

namespace walker {

void
registerWrightFisher( DiffEqFactory& f, std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register Wright-Fisher SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for SDE
  using WrightFisherPolicies =
    tk::cartesian_product< InitPolicies, WrightFisherCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< WrightFisherPolicies >(
    registerDiffEq< WrightFisher >( f, ctr::DiffEqType::WRIGHTFISHER, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoWrightFisher( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the Wright-Fisher SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::WRIGHTFISHER ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::WRIGHTFISHER ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::wrightfisher >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::wrightfisher >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::wrightfisher, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff omega [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::wrightfisher, tag::omega >().at(c) ) );
  spikes( nfo, g_inputdeck.get< tag::param, tag::wrightfisher, tag::init,
                                tag::spike >().at(c)
  );
  betapdfs( nfo, g_inputdeck.get< tag::param, tag::wrightfisher, tag::init,
                                  tag::betapdf >().at(c) );

  return nfo;
}

}  // walker::
