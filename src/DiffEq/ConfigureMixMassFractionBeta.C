// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureMixMassFractionBeta.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
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

#include "Tags.h"
#include "CartesianProduct.h"
#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"
#include "Walker/Options/InitPolicy.h"

#include "ConfigureMixMassFractionBeta.h"
#include "MixMassFractionBeta.h"
#include "MixMassFractionBetaCoeffPolicy.h"

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
infoMixMassFractionBeta( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt )
// *****************************************************************************
//  Return information on the mix mass fraction beta SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::MIXMASSFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back(
    ctr::DiffEq().name( ctr::DiffEqType::MIXMASSFRACBETA ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::mixmassfracbeta >(c) ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::mixmassfracbeta >()[c];
  nfo.emplace_back( "number of components",
    std::to_string( ncomp ) + " (" + std::to_string(ncomp/4) + "*4) " );

  // Optional coupled
  const auto& coupled_velocity =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::velocity >();
  const auto& coupled_velocity_id =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::velocity_id >();
  Assert( coupled_velocity.size() == coupled_velocity_id.size(),
          "Size mismatch" );

  if (coupled_velocity.size() > c) {
    nfo.emplace_back( "coupled velocity depvar", std::string( 1,
      coupled_velocity[c] ) );
    nfo.emplace_back( "coupled velocity depvar offset", std::to_string(
      coupled_velocity_id[c] ) );
  }

  // Optional coupled
  const auto& coupled_dissipation =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::dissipation >();
  const auto& coupled_dissipation_id =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::dissipation_id >();
  Assert( coupled_dissipation.size() == coupled_dissipation_id.size(),
          "Size mismatch" );

  if (coupled_dissipation.size() > c) {
    nfo.emplace_back( "coupled dissipation depvar", std::string( 1,
      coupled_dissipation[c] ) );
    nfo.emplace_back( "coupled dissipation depvar ofs", std::to_string(
    coupled_dissipation_id[c] ) );
  }

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::initpolicy >()[c] )
  );
  auto cp =
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::coeffpolicy >()[c];
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name( cp ) );
  if (cp == ctr::CoeffPolicyType::HYDROTIMESCALE) {
    nfo.emplace_back(
      "inverse hydro time scales [" + std::to_string( ncomp/4 ) + "]",
      options( ctr::HydroTimeScales(),
               g_inputdeck.get< tag::param,
                                tag::mixmassfracbeta,
                                tag::hydrotimescales >().at(c) ) );
    nfo.emplace_back(
      "production/dissipation [" + std::to_string( ncomp/4 ) + "]",
      options( ctr::HydroProductions(),
               g_inputdeck.get< tag::param,
                                tag::mixmassfracbeta,
                                tag::hydroproductions >().at(c) ) );
  }

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b' [" + std::to_string( ncomp/4 ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::bprime >().at(c) )
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp/4 ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::S >().at(c) )
  );
  nfo.emplace_back(
    "coeff kappa' [" + std::to_string( ncomp/4 ) + "]",
    parameters( g_inputdeck.get< tag::param,
                                 tag::mixmassfracbeta,
                                 tag::kappaprime >().at(c) ) );
  nfo.emplace_back(
    "coeff rho2 [" + std::to_string( ncomp/4 ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::rho2 >().at(c) ) );
  nfo.emplace_back(
    "coeff r [" + std::to_string( ncomp/4 ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::r >().at(c) )
  );
  spikes( nfo,
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::spike >().at(c) );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::mixmassfracbeta, tag::betapdf >().at(c) );

  return nfo;
}

}  // walker::
