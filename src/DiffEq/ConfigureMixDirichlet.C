// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureMixDirichlet.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the MixDirichlet SDE
  \details   Register and compile configuration on the MixDirichlet SDE.
*/
// *****************************************************************************

#include <set>
#include <map>
#include <vector>
#include <string>
#include <utility>

#include <brigand/algorithms/for_each.hpp>

#include "Tags.h"
#include "CartesianProduct.h"
#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"
#include "Walker/Options/InitPolicy.h"

#include "ConfigureMixDirichlet.h"
#include "MixDirichlet.h"
#include "MixDirichletCoeffPolicy.h"

namespace walker {

void
registerMixDirichlet( DiffEqFactory& f, std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register MixDirichlet SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for SDE
  using MixDirPolicies =
    tk::cartesian_product< InitPolicies, MixDirichletCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< MixDirPolicies >(
    registerDiffEq< MixDirichlet >( f, ctr::DiffEqType::MIXDIRICHLET, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoMixDirichlet( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt )
// *****************************************************************************
//  Return information on the MixDirichlet SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  using eq = tag::mixdirichlet;

  auto c = ++cnt[ ctr::DiffEqType::MIXDIRICHLET ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::MIXDIRICHLET ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< eq >(c) ) );
  const auto ncomp = g_inputdeck.get< tag::component >().get< eq >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, eq, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, eq, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, eq, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, eq, tag::rng >()[c] ) );

  auto K = ncomp - MIXDIR_NUMDERIVED;
  auto N = K + 1;

  nfo.emplace_back( "coeff b [" + std::to_string(K) + "]",
    parameters( g_inputdeck.get< tag::param, eq, tag::b >().at(c) )
  );
  nfo.emplace_back( "coeff S [" + std::to_string(K) + "]",
    parameters( g_inputdeck.get< tag::param, eq, tag::S >().at(c) )
  );
  nfo.emplace_back( "coeff kappaprime [" + std::to_string(K) + "]",
    parameters( g_inputdeck.get< tag::param, eq, tag::kappaprime >().at(c) )
  );

  const auto& rho = g_inputdeck.get< tag::param, eq, tag::rho >();
  if (!rho.empty()) {
    nfo.emplace_back( "coeff rho [" + std::to_string(N) + "]",
                      parameters( rho.at(c) ) );
    nfo.emplace_back( "coeff r [" + std::to_string(K) + "]",
                      parameters( MixDir_r(rho[c]) ) );
  }

  return nfo;
}

}  // walker::
