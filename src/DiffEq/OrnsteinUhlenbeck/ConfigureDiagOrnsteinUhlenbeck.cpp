// *****************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeck/ConfigureDiagOrnsteinUhlenbeck.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the diagonal
             Ornstein-Uhlenbeck SDE
  \details   Register and compile configuration on the diagonal
             Ornstein-Uhlenbeck SDE.
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

#include "ConfigureDiagOrnsteinUhlenbeck.hpp"
#include "DiagOrnsteinUhlenbeck.hpp"
#include "DiagOrnsteinUhlenbeckCoeffPolicy.hpp"

namespace walker {

void
registerDiagOrnsteinUhlenbeck( DiffEqFactory& f,
                               std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register diagonal Ornstein-Uhlenbeck SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for SDE
  using DiagOrnsteinUhlenbeckPolicies =
    tk::cartesian_product< InitPolicies, DiagOrnsteinUhlenbeckCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< DiagOrnsteinUhlenbeckPolicies >(
    registerDiffEq< DiagOrnsteinUhlenbeck >( f, ctr::DiffEqType::DIAG_OU, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoDiagOrnsteinUhlenbeck( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the diagonal Ornstein-Uhlenbeck SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::DIAG_OU ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DIAG_OU ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::diagou >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::diagou >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::diagou, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::diagou, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff sigmasq [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::diagou, tag::sigmasq >().at(c) ) );
  nfo.emplace_back( "coeff theta [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::diagou, tag::theta >().at(c) )
  );
  nfo.emplace_back( "coeff mu [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::diagou, tag::mu >().at(c) ) );
  spikes( nfo,
          g_inputdeck.get< tag::param, tag::diagou, tag::spike >().at(c) );
  betapdfs( nfo,
            g_inputdeck.get< tag::param, tag::diagou, tag::betapdf >().at(c) );

  return nfo;
}

}  // walker::
