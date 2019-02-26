// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureSkewNormal.C
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Register and compile configuration on the skew-normal SDE
  \details   Register and compile configuration on the skew-normal SDE.
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

#include "ConfigureSkewNormal.h"
#include "SkewNormal.h"
#include "SkewNormalCoeffPolicy.h"

namespace walker {

void
registerSkewNormal( DiffEqFactory& f, std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register skew-normal SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for SDE
  using SkewNormalPolicies =
    tk::cartesian_product< InitPolicies, SkewNormalCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< SkewNormalPolicies >(
    registerDiffEq< SkewNormal >( f, ctr::DiffEqType::SKEWNORMAL, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoSkewNormal( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt )
// *****************************************************************************
//  Return information on the skew-normal SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::SKEWNORMAL ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::SKEWNORMAL ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::skewnormal >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::skewnormal >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::skewnormal, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::skewnormal, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff T [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::skewnormal, tag::timescale >().at(c) )
  );
  nfo.emplace_back(
    "coeff sigmasq [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::skewnormal, tag::sigmasq >().at(c) ) );
  nfo.emplace_back(
    "coeff lambda [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::skewnormal, tag::lambda >().at(c) ) );
  spikes( nfo,
          g_inputdeck.get< tag::param, tag::skewnormal, tag::spike >().at(c) );
  betapdfs(
    nfo,
    g_inputdeck.get< tag::param, tag::skewnormal, tag::betapdf >().at(c) );

  return nfo;
}

}  // walker::
