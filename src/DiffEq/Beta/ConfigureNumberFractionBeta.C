// *****************************************************************************
/*!
  \file      src/DiffEq/Beta/ConfigureNumberFractionBeta.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the number fraction beta SDE
  \details   Register and compile configuration on the number fraction beta SDE.
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

#include "ConfigureNumberFractionBeta.h"
#include "NumberFractionBeta.h"
#include "NumberFractionBetaCoeffPolicy.h"

namespace walker {

void
registerNumberFractionBeta( DiffEqFactory& f, std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register number fraction beta SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for SDE
  using NumberFractionBetaPolicies =
    tk::cartesian_product< InitPolicies, NumberFractionBetaCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< NumberFractionBetaPolicies >(
    registerDiffEq< NumberFractionBeta >
                  ( f, ctr::DiffEqType::NUMFRACBETA, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoNumberFractionBeta( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt )
// *****************************************************************************
//  Return information on the number fraction beta SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::NUMFRACBETA ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::NUMFRACBETA ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::numfracbeta >(c) ) );
  auto ncomp =
    g_inputdeck.get< tag::component >().get< tag::numfracbeta >()[c] / 3;
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::coeffpolicy >()[c] ) );

  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::numfracbeta, tag::b >().at(c) )
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::numfracbeta, tag::S >().at(c) )
  );
  nfo.emplace_back(
    "coeff kappa [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::numfracbeta, tag::kappa >().at(c) ) );
  nfo.emplace_back(
    "coeff rho2 [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::numfracbeta, tag::rho2 >().at(c) ) );
  nfo.emplace_back(
    "coeff rcomma [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::numfracbeta, tag::rcomma >().at(c) ) );
  spikes( nfo,
          g_inputdeck.get< tag::param, tag::numfracbeta, tag::spike >().at(c) );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::numfracbeta, tag::betapdf >().at(c) );

  return nfo;
}

}  // walker::
