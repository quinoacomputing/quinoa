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
#include "InfoMesh.hpp"

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
  using tag::param;
  using tk::parameters;
  using eq = tag::transport;

  auto c = ++cnt[ ctr::PDEType::TRANSPORT ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::TRANSPORT ), "" );

  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< param, eq, tag::depvar >()[c] ) );

  infoMesh< eq >( c, nfo );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< param, eq, tag::problem >()[c] ) );

  auto intsharp = g_inputdeck.get< tag::param, eq, tag::intsharp >()[c];
  nfo.emplace_back( "interface sharpening", std::to_string( intsharp ) );

  if (intsharp)
  {
    auto intsharp_param =
      g_inputdeck.get< tag::param, eq, tag::intsharp_param >()[c];
    nfo.emplace_back( "interface sharpening parameter",
                      std::to_string( intsharp_param ) );
  }

  auto ncomp = g_inputdeck.get< tag::component >().get< eq >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  const auto& diff = g_inputdeck.get< param, eq, tag::diffusivity >();
  if (diff.size() > c)
    nfo.emplace_back( "coeff diffusivity [" + std::to_string( ncomp ) + "]",
                       parameters( diff[c] ) );

  const auto& u0 = g_inputdeck.get< param, eq, tag::u0 >();
  if (u0.size() > c)
    nfo.emplace_back( "coeff u0 [" + std::to_string( ncomp ) + "]",
                       parameters( u0[c] ) );

  const auto& lambda = g_inputdeck.get< param, eq, tag::lambda >();
  if (lambda.size() > c)
    nfo.emplace_back( "coeff lambda [" + std::to_string( ncomp ) + "]",
      parameters( lambda[c] ) );

  const auto& bcdir = g_inputdeck.get< param, eq, tag::bc, tag::bcdir >();
  if (bcdir.size() > c)
    nfo.emplace_back( "Dirichlet boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcdir[c] ) );

  const auto& bcsym = g_inputdeck.get< param, eq, tag::bc, tag::bcsym >();
  if (bcsym.size() > c)
    nfo.emplace_back( "Symmetry boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcsym[c] ) );

  const auto& bcinlet =
    g_inputdeck.get< param, eq, tag::bc, tag::bcinlet >();
  if (bcinlet.size() > c)
    nfo.emplace_back( "Inlet boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcinlet[c] ) );

  const auto& bcoutlet =
    g_inputdeck.get< param, eq, tag::bc, tag::bcoutlet >();
  if (bcoutlet.size() > c)
    nfo.emplace_back( "Outlet boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcoutlet[c] ) );

  const auto& bcextrapolate =
    g_inputdeck.get< param, eq, tag::bc, tag::bcextrapolate >();
  if (bcextrapolate.size() > c)
    nfo.emplace_back( "Symmetry boundary [" + std::to_string( ncomp ) + "]",
      parameters( bcextrapolate[c] ) );

  return nfo;
}

void
assignTransportGetVars( const std::string& name, tk::GetVarFn& f )
// *****************************************************************************
// Assign functions that compute physics variables from the numerical solution
// for Transport
//! \param[in] name Name of variable whose tk::GetVarFn is to be assigned
//! \param[in,out] f Function assigned
// *****************************************************************************
{
  using namespace kw;
  using namespace transport;

  assign< outvar_material_indicator >( name, matIndicatorOutVar, f );
}

}  // inciter::
