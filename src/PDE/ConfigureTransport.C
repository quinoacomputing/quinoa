// *****************************************************************************
/*!
  \file      src/PDE/ConfigureTransport.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Register and compile configuration on the transport PDE
  \details   Register and compile configuration on the transport PDE.
*/
// *****************************************************************************

#include <set>
#include <map>
#include <vector>
#include <string>

#include <brigand/algorithms/for_each.hpp>

#include "Tags.h"
#include "CartesianProduct.h"
#include "PDEFactory.h"
#include "Inciter/Options/PDE.h"

#include "ConfigureTransport.h"
#include "Transport/Physics/CG.h"
#include "Transport/Physics/DG.h"
#include "Transport/CGTransport.h"
#include "Transport/DGTransport.h"
#include "Transport/Problem.h"

namespace inciter {

void
registerTransport( CGFactory& cf,
                   DGFactory& df,
                   std::set< ctr::PDEType >& cgt,
                   std::set< ctr::PDEType >& dgt )
// *****************************************************************************
// Register transport PDE into PDE factory
//! \param[in,out] cf Continuous Galerkin PDE factory to register to
//! \param[in,out] df Discontinuous Galerkin PDE factory to register to
//! \param[in,out] cgt Counters for equation types registered into CG factory
//! \param[in,out] dgt Counters for equation types registered into DG factory
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies
  using CGTransportPolicies =
    tk::cartesian_product< cg::TransportPhysics, TransportProblems >;
  // Register PDEs for all combinations of policies
  brigand::for_each< CGTransportPolicies >(
    registerCG< cg::Transport >( cf, cgt, ctr::PDEType::TRANSPORT ) );

  // Construct vector of vectors for all possible policies
  using DGTransportPolicies =
    tk::cartesian_product< dg::TransportPhysics, TransportProblems >;
  // Register PDEs for all combinations of policies
  brigand::for_each< DGTransportPolicies >(
    registerDG< dg::Transport >( df, dgt, ctr::PDEType::TRANSPORT ) );
}

std::vector< std::pair< std::string, std::string > >
infoTransport( std::map< ctr::PDEType, tk::ctr::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the transport PDE
//! \param[inout] cnt std::map of counters for all PDE types
//! \return vector of string pairs describing the PDE configuration
// *****************************************************************************
{
  auto c = ++cnt[ ctr::PDEType::TRANSPORT ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::TRANSPORT ), "" );

  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, tag::transport, tag::depvar >()[c] ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< tag::param, tag::transport, tag::problem >()[c] ) );

  nfo.emplace_back( "start offset in unknowns array", std::to_string(
    g_inputdeck.get< tag::component >().offset< tag::transport >(c) ) );

  auto ncomp = g_inputdeck.get< tag::component >().get< tag::transport >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  const auto& diff =
     g_inputdeck.get< tag::param, tag::transport, tag::diffusivity >();
  if (diff.size() > c)
    nfo.emplace_back( "coeff diffusivity [" + std::to_string( ncomp ) + "]",
                       parameters( diff[c] ) );

  const auto& u0 = g_inputdeck.get< tag::param, tag::transport, tag::u0 >();
  if (u0.size() > c)
    nfo.emplace_back( "coeff u0 [" + std::to_string( ncomp ) + "]",
                       parameters( u0[c] ) );

  const auto& lambda =
    g_inputdeck.get< tag::param, tag::transport, tag::lambda >();
  if (lambda.size() > c)
    nfo.emplace_back( "coeff lambda [" + std::to_string( ncomp ) + "]",
      parameters( lambda[c] ) );

  const auto& bcdir =
    g_inputdeck.get< tag::param, tag::transport, tag::bcdir >();
  if (bcdir.size() > c)
    nfo.emplace_back( "Dirichlet boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcdir[c] ) );

  const auto& bcsym =
    g_inputdeck.get< tag::param, tag::transport, tag::bcsym >();
  if (bcsym.size() > c)
    nfo.emplace_back( "Symmetry boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcsym[c] ) );

  const auto& bcinlet =
    g_inputdeck.get< tag::param, tag::transport, tag::bcinlet >();
  if (bcinlet.size() > c)
    nfo.emplace_back( "Inlet boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcinlet[c] ) );

  const auto& bcoutlet =
    g_inputdeck.get< tag::param, tag::transport, tag::bcoutlet >();
  if (bcoutlet.size() > c)
    nfo.emplace_back( "Outlet boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcoutlet[c] ) );

  const auto& bcextrapolate =
    g_inputdeck.get< tag::param, tag::transport, tag::bcextrapolate >();
  if (bcextrapolate.size() > c)
    nfo.emplace_back( "Symmetry boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcextrapolate[c] ) );

  return nfo;
}

}  // inciter::
