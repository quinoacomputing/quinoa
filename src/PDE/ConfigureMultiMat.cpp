// *****************************************************************************
/*!
  \file      src/PDE/ConfigureMultiMat.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration for multi-material compressible
     flow PDE
  \details   Register and compile configuration for compressible multi-material
     flow PDE.
*/
// *****************************************************************************

#include <set>
#include <map>
#include <vector>
#include <string>

#include <brigand/algorithms/for_each.hpp>

#include "Tags.hpp"
#include "CartesianProduct.hpp"
#include "PDEFactory.hpp"
#include "Inciter/Options/PDE.hpp"

#include "ConfigureMultiMat.hpp"
#include "MultiMat/Physics/DG.hpp"
#include "MultiMat/DGMultiMat.hpp"
#include "MultiMat/Problem.hpp"

namespace inciter {

void
registerMultiMat( DGFactory& df, std::set< ctr::PDEType >& dgt )
// *****************************************************************************
// Register multi-material compressible flow PDE into PDE factory
//! \param[in,out] df Discontinuous Galerkin PDE factory to register to
//! \param[in,out] dgt Counters for equation types registered into DG factory
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies
  using DGMultiMatPolicies =
    tk::cartesian_product< dg::MultiMatPhysics, MultiMatProblems >;
  // Register PDEs for all combinations of policies
  brigand::for_each< DGMultiMatPolicies >(
    registerDG< dg::MultiMat >( df, dgt, ctr::PDEType::MULTIMAT ) );
}

std::vector< std::pair< std::string, std::string > >
infoMultiMat( std::map< ctr::PDEType, tk::ctr::ncomp_type >& cnt )
// *****************************************************************************
//  Return information on the compressible flow system of PDEs
//! \param[inout] cnt std::map of counters for all PDE types
//! \return vector of string pairs describing the PDE configuration
// *****************************************************************************
{
  using eq = tag::multimat;

  auto c = ++cnt[ ctr::PDEType::MULTIMAT ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::MULTIMAT ), "" );

  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, eq, tag::depvar >()[c] ) );

  nfo.emplace_back( "physics", ctr::Physics().name(
    g_inputdeck.get< tag::param, eq, tag::physics >()[c] ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< tag::param, eq, tag::problem >()[c] ) );

  auto nmat = g_inputdeck.get< tag::param, eq, tag::nmat >()[c];
  nfo.emplace_back( "number of materials", std::to_string( nmat ) );

  auto ncomp = g_inputdeck.get< tag::component >().get< eq >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  nfo.emplace_back( "start offset in unknowns array", std::to_string(
    g_inputdeck.get< tag::component >().offset< eq >(c) ) );

  nfo.emplace_back( "material id", parameters(
    g_inputdeck.get< tag::param, eq, tag::id >() ) );

  nfo.emplace_back( "ratio of specific heats", parameters(
    g_inputdeck.get< tag::param, eq, tag::gamma >() ) );

  const auto& alpha = g_inputdeck.get< tag::param, eq, tag::alpha >();
  if (!alpha.empty()) nfo.emplace_back( "coeff alpha", parameters( alpha ) );

  const auto& beta = g_inputdeck.get< tag::param, eq, tag::beta >();
  if (!beta.empty())
    nfo.emplace_back( "coeff beta", parameters( beta ) );

  const auto& p0 = g_inputdeck.get< tag::param, eq, tag::p0 >();
  if (!p0.empty())
    nfo.emplace_back( "coeff p0", parameters( p0 ) );

  return nfo;
}

}  // inciter::
