// *****************************************************************************
/*!
  \file      src/DiffEq/Beta/ConfigureMixMassFractionBeta.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the mix mass fraction beta
             SDE
  \details   Register and compile configuration on the mix mass fraction beta
             SDE.
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

#include "ConfigureMixMassFractionBeta.hpp"
#include "MixMassFractionBeta.hpp"
#include "MixMassFractionBetaCoeffPolicy.hpp"

#include "CoupledEq.hpp"

namespace walker {

void
registerMixMassFractionBeta( DiffEqFactory& f, std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register mix mass fraction beta SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for SDE
  using MixMassFracBetaPolicies =
    tk::cartesian_product< InitPolicies, MixMassFracBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< MixMassFracBetaPolicies >(
    registerDiffEq< MixMassFractionBeta >
                  ( f, ctr::DiffEqType::MIXMASSFRACBETA, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoMixMassFractionBeta( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the mix mass fraction beta SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  using eq = tag::mixmassfracbeta;

  auto c = ++cnt[ ctr::DiffEqType::MIXMASSFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back(
    ctr::DiffEq().name( ctr::DiffEqType::MIXMASSFRACBETA ), "" );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, eq, tag::depvar >()[c] ) );
  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< eq >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component, eq >()[c];

  auto numderived =
    MixMassFractionBeta<InitZero,MixMassFracBetaCoeffInstVel>::NUMDERIVED;
  nfo.emplace_back( "number of components", std::to_string(ncomp) + " (=" +
                    std::to_string(ncomp/(numderived+1)) + '*' +
                    std::to_string(numderived+1) + ") " );

  coupledInfo< eq, tag::velocity, tag::velocity_id >
             ( c, "velocity", nfo );
  coupledInfo< eq, tag::dissipation, tag::dissipation_id >
             ( c, "dissipation", nfo );

  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, eq, tag::initpolicy >()[c] )
  );
  auto cp = g_inputdeck.get< tag::param, eq, tag::coeffpolicy >()[c];
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name( cp ) );
  if (cp == ctr::CoeffPolicyType::HYDROTIMESCALE) {
    nfo.emplace_back(
      "inverse hydro time scales [" + std::to_string( ncomp/4 ) + "]",
      options( ctr::HydroTimeScales(),
               g_inputdeck.get< tag::param,
                                eq,
                                tag::hydrotimescales >().at(c) ) );
    nfo.emplace_back(
      "production/dissipation [" + std::to_string( ncomp/4 ) + "]",
      options( ctr::HydroProductions(),
               g_inputdeck.get< tag::param,
                                eq,
                                tag::hydroproductions >().at(c) ) );
  }

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, eq, tag::rng >()[c] ) );

  const auto& mg = g_inputdeck.get< tag::param, eq, tag::mean_gradient >();
  if (mg.size() > c) {
    const auto& mean_gradient = mg[c];
    Assert( mean_gradient.size() == 3, "Mean gradient vector size must be 3" );
    nfo.emplace_back( "imposed mean gradient", parameters(mean_gradient) );
  }

  nfo.emplace_back(
    "coeff b' [" + std::to_string( ncomp/4 ) + "]",
    parameters( g_inputdeck.get< tag::param, eq, tag::bprime >().at(c) )
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp/4 ) + "]",
    parameters( g_inputdeck.get< tag::param, eq, tag::S >().at(c) )
  );
  nfo.emplace_back(
    "coeff kappa' [" + std::to_string( ncomp/4 ) + "]",
    parameters( g_inputdeck.get< tag::param,
                                 eq,
                                 tag::kappaprime >().at(c) ) );
  nfo.emplace_back(
    "coeff rho2 [" + std::to_string( ncomp/4 ) + "]",
    parameters( g_inputdeck.get< tag::param, eq, tag::rho2 >().at(c) ) );
  nfo.emplace_back(
    "coeff r [" + std::to_string( ncomp/4 ) + "]",
    parameters( g_inputdeck.get< tag::param, eq, tag::r >().at(c) )
  );

  spikes( nfo,
    g_inputdeck.get< tag::param, eq, tag::init, tag::spike >().at(c) );
  betapdfs( nfo,
    g_inputdeck.get< tag::param, eq, tag::init, tag::betapdf >().at(c) );

  return nfo;
}

}  // walker::
