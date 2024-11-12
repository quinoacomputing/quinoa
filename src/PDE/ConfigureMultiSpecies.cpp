// *****************************************************************************
/*!
  \file      src/PDE/ConfigureMultiSpecies.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration for multi-species compressible
     flow PDE
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
#include "ConfigureMultiSpecies.hpp"
#include "MultiSpecies/Physics/DG.hpp"
#include "MultiSpecies/DGMultiSpecies.hpp"
#include "MultiSpecies/Problem.hpp"
#include "Inciter/Options/Material.hpp"

namespace inciter {

void
registerMultiSpecies( DGFactory& df, FVFactory& /*ff*/,
  std::set< ctr::PDEType >& /*fvt*/, std::set< ctr::PDEType >& dgt )
// *****************************************************************************
// Register multi-material compressible flow PDE into PDE factory
//! \param[in,out] df Discontinuous Galerkin PDE factory to register to
// //! \param[in,out] ff Finite volume PDE factory to register to
//! \param[in,out] dgt Counters for equation types registered into DG factory
// //! \param[in,out] fvt Counters for equation types registered into FV factory
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies
  using DGMultiSpeciesPolicies =
    tk::cartesian_product< dg::MultiSpeciesPhysics, MultiSpeciesProblems >;
  // Register PDEs for all combinations of policies
  brigand::for_each< DGMultiSpeciesPolicies >(
    registerDG< dg::MultiSpecies >( df, dgt, ctr::PDEType::MULTISPECIES ) );
}

std::vector< std::pair< std::string, std::string > >
infoMultiSpecies( std::map< ctr::PDEType, tk::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the compressible flow system of PDEs
//! \param[inout] cnt std::map of counters for all PDE types
//! \return vector of string pairs describing the PDE configuration
// *****************************************************************************
{
  using eq = tag::multispecies;
  using tk::parameter;
  using tk::parameters;

  auto c = ++cnt[ ctr::PDEType::MULTISPECIES ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::MULTISPECIES ), "" );

  nfo.emplace_back( "physics", ctr::Physics().name(
    g_inputdeck.get< eq, tag::physics >() ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< eq, tag::problem >() ) );

  nfo.emplace_back( "flux", ctr::Flux().name(
    g_inputdeck.get< tag::flux >() ) );

  auto nspec = g_inputdeck.get< eq, tag::nspec >();
  nfo.emplace_back( "number of species", std::to_string( nspec ) );

  auto viscous = g_inputdeck.get< eq, tag::viscous >();
  nfo.emplace_back( "viscosity", std::to_string( viscous ) );

  auto ncomp = g_inputdeck.get< tag::ncomp >();
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  // ICs and IC-boxes

  const auto& ic = g_inputdeck.get< tag::ic >();

  const auto& bgmatidic = ic.get< tag::materialid >();
  nfo.emplace_back( "IC background material id", parameter( bgmatidic ) );

  const auto& icbox = ic.get< tag::box >();
  if (!icbox.empty()) {
    std::size_t bcnt = 0;
    for (const auto& b : icbox) {   // for all boxes configured for this eq
      std::vector< tk::real > box
        { b.get< tag::xmin >(), b.get< tag::xmax >(),
          b.get< tag::ymin >(), b.get< tag::ymax >(),
          b.get< tag::zmin >(), b.get< tag::zmax >() };

      std::string boxname = "IC box " + parameter(bcnt);
      nfo.emplace_back( boxname, parameters( box ) );

      nfo.emplace_back( boxname + " orientation",
        parameters(b.get< tag::orientation >()) );

      nfo.emplace_back( boxname + " mass fractions",
        parameters(b.get< tag::mass_fractions >()) );

      const auto& initiate = b.get< tag::initiate >();
      auto opt = ctr::Initiate();
      nfo.emplace_back( boxname + ' ' + opt.group(), opt.name(initiate) );

      ++bcnt;
    }
  }

  const auto& icblock = ic.get< tag::meshblock >();
  for (const auto& b : icblock) {   // for all blocks configured for eq
    std::string blockname = "IC mesh block " +
      parameter(b.get< tag::blockid >());

    nfo.emplace_back( blockname + " material id",
                      parameter( b.get< tag::materialid >() ) );
    const auto& initiate = b.get< tag::initiate >();
    auto opt = ctr::Initiate();
    nfo.emplace_back( blockname + ' ' + opt.group(), opt.name(initiate) );
  }

  return nfo;
}

}  // inciter::
