// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureDirichlet.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Register and compile configuration on the Dirichlet SDE
  \details   Register and compile configuration on the Dirichlet SDE.
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

#include "ConfigureDirichlet.h"
#include "Dirichlet.h"
#include "DirichletCoeffPolicy.h"

namespace walker {

void
registerDirichlet( DiffEqFactory& f, std::set< ctr::DiffEqType >& t )
// *****************************************************************************
// Register Dirichlet SDE into DiffEq factory
//! \param[in,out] f Differential equation factory to register to
//! \param[in,out] t Counters for equation types registered
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies for SDE
  using DirPolicies =
    tk::cartesian_product< InitPolicies, DirichletCoeffPolicies >;
  // Register SDE for all combinations of policies
  brigand::for_each< DirPolicies >(
    registerDiffEq< Dirichlet >( f, ctr::DiffEqType::DIRICHLET, t ) );
}

std::vector< std::pair< std::string, std::string > >
infoDirichlet( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt )
// *****************************************************************************
//  Return information on the Dirichlet SDE
//! \param[inout] cnt std::map of counters for all differential equation types
//! \return vector of string pairs describing the SDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::DiffEqType::DIRICHLET ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::DiffEq().name( ctr::DiffEqType::DIRICHLET ), "" );

  nfo.emplace_back( "start offset in particle array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::dirichlet >(c) ) );
  auto ncomp = g_inputdeck.get< tag::component >().get< tag::dirichlet >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "kind", "stochastic" );
  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::dirichlet, tag::depvar >()[c] ) );
  nfo.emplace_back( "initialization policy", ctr::InitPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::initpolicy >()[c] ) );
  nfo.emplace_back( "coefficients policy", ctr::CoeffPolicy().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::coeffpolicy >()[c] ) );
  nfo.emplace_back( "random number generator", tk::ctr::RNG().name(
    g_inputdeck.get< tag::param, tag::dirichlet, tag::rng >()[c] ) );
  nfo.emplace_back(
    "coeff b [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::dirichlet, tag::b >().at(c) )
  );
  nfo.emplace_back(
    "coeff S [" + std::to_string( ncomp ) + "]",
    parameters( g_inputdeck.get< tag::param, tag::dirichlet, tag::S >().at(c) )
  );
  nfo.emplace_back(
    "coeff kappa [" + std::to_string( ncomp ) + "]",
    parameters(
      g_inputdeck.get< tag::param, tag::dirichlet, tag::kappa >().at(c) )
  );

  return nfo;
}

}  // walker::
