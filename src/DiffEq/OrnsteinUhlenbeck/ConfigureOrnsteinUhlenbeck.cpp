// *****************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeck/ConfigureOrnsteinUhlenbeck.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the Ornstein-Uhlenbeck SDE
  \details   Register and compile configuration on the Ornstein-Uhlenbeck SDE.
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

#include "ConfigureOrnsteinUhlenbeck.hpp"
#include "OrnsteinUhlenbeck.hpp"
#include "OrnsteinUhlenbeckCoeffPolicy.hpp"

namespace walker {

void
registerOrnsteinUhlenbeck( DiffEqFactory& f, std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register Ornstein-Uhlenbeck SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for SDE
  using OrnsteinUhlenbeckPolicies =
    tk::cartesian_product< InitPolicies, OrnsteinUhlenbeckCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< OrnsteinUhlenbeckPolicies >(
    registerDiffEq< OrnsteinUhlenbeck >( f, ctr::DiffEqType::OU, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoOrnsteinUhlenbeck( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt )
// *****************************************************************************
//  Return information on the Ornstein-Uhlenbeck SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::OU ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::OU ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::ou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::ou >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::ou, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::ou, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::ou, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff sigmasq [" + std::to_string( ncomp*(ncomp+1)/2 ) + ", upper tri]",
    parameters( g_inputdeck.get< tag::param, tag::ou, tag::sigmasq >().at(c) )
  );
  nfo.emplace_back( "coeff theta [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::ou, tag::theta >().at(c) ) );
  nfo.emplace_back( "coeff mu [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::ou, tag::mu >().at(c) ) );
  spikes( nfo, g_inputdeck.get< tag::param, tag::ou, tag::init,
                                tag::spike >().at(c) );
  betapdfs( nfo, g_inputdeck.get< tag::param, tag::ou, tag::init,
                                  tag::betapdf >().at(c) );

  return nfo;
}

}  // walker::
