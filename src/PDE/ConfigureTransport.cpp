// *****************************************************************************
/*!
  \file      src/PDE/ConfigureTransport.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the transport PDE
  \details   Register and compile configuration on the transport PDE.
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
#include "ContainerUtil.hpp"
#include "ConfigureTransport.hpp"
#include "Transport/Physics/CG.hpp"
#include "Transport/Physics/DG.hpp"
#include "Transport/CGTransport.hpp"
#include "Transport/DGTransport.hpp"
#include "Transport/Problem.hpp"

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
infoTransport( std::map< ctr::PDEType, tk::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the transport PDE
//! \param[inout] cnt std::map of counters for all PDE types
//! \return vector of string pairs describing the PDE configuration
// *****************************************************************************
{
  using tk::parameters;
  using eq = tag::transport;

  auto c = ++cnt[ ctr::PDEType::TRANSPORT ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::TRANSPORT ), "" );

  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::depvar >()[c] ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< eq, tag::problem >() ) );

  auto intsharp = g_inputdeck.get< eq, tag::intsharp >();
  nfo.emplace_back( "interface sharpening", std::to_string( intsharp ) );

  auto ncomp = g_inputdeck.get< tag::ncomp >();
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  const auto& bc = g_inputdeck.get< tag::bc >();
  for (const auto& ib : bc) {
    const auto& bcdir = ib.get< tag::dirichlet >();
    if (!bcdir.empty())
      nfo.emplace_back( "Dirichlet boundary [" + std::to_string( ncomp ) + "]",
        parameters( bcdir ) );

    const auto& bcsym = ib.get< tag::symmetry >();
    if (!bcsym.empty())
      nfo.emplace_back( "Symmetry boundary [" + std::to_string( ncomp ) + "]",
        parameters( bcsym ) );

    const auto& bcinlet =
      ib.get< tag::inlet >();
    if (!bcinlet.empty())
      nfo.emplace_back( "Inlet boundary [" + std::to_string( ncomp ) + "]",
        parameters( bcinlet ) );

    const auto& bcoutlet =
      ib.get< tag::outlet >();
    if (!bcoutlet.empty())
      nfo.emplace_back( "Outlet boundary [" + std::to_string( ncomp ) + "]",
        parameters( bcoutlet ) );

    const auto& bcextrapolate =
      ib.get< tag::extrapolate >();
    if (!bcextrapolate.empty())
      nfo.emplace_back( "Symmetry boundary [" + std::to_string( ncomp ) + "]",
        parameters( bcextrapolate ) );
  }

  return nfo;
}

}  // inciter::
