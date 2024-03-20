// *****************************************************************************
/*!
  \file      src/PDE/ConfigureCompFlow.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration for compressible flow PDE
  \details   Register and compile configuration for compressible flow PDE.
*/
// *****************************************************************************

#include <set>
#include <map>
#include <vector>
#include <string>
#include <limits>

#include <brigand/algorithms/for_each.hpp>

#include "Tags.hpp"
#include "CartesianProduct.hpp"
#include "PDEFactory.hpp"
#include "Inciter/Options/PDE.hpp"
#include "ContainerUtil.hpp"
#include "ConfigureCompFlow.hpp"
#include "CompFlow/Physics/CG.hpp"
#include "CompFlow/Physics/DG.hpp"
#include "CompFlow/CGCompFlow.hpp"
#include "CompFlow/DGCompFlow.hpp"
#include "CompFlow/Problem.hpp"

namespace inciter {

void
registerCompFlow( CGFactory& cf,
                  DGFactory& df,
                  std::set< ctr::PDEType >& cgt,
                  std::set< ctr::PDEType >& dgt )
// *****************************************************************************
// Register compressible flow PDE into PDE factory
//! \param[in,out] cf Continuous Galerkin PDE factory to register to
//! \param[in,out] df Discontinuous Galerkin PDE factory to register to
//! \param[in,out] cgt Counters for equation types registered into CG factory
//! \param[in,out] dgt Counters for equation types registered into DG factory
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies
  using CGCompFlowPolicies =
    tk::cartesian_product< cg::CompFlowPhysics, CompFlowProblems >;
  // Register PDEs for all combinations of policies
  brigand::for_each< CGCompFlowPolicies >(
    registerCG< cg::CompFlow >( cf, cgt, ctr::PDEType::COMPFLOW ) );

  // Construct vector of vectors for all possible policies
  using DGCompFlowPolicies =
    tk::cartesian_product< dg::CompFlowPhysics, CompFlowProblems >;
  // Register PDEs for all combinations of policies
  brigand::for_each< DGCompFlowPolicies >(
    registerDG< dg::CompFlow >( df, dgt, ctr::PDEType::COMPFLOW ) );
}

std::vector< std::pair< std::string, std::string > >
infoCompFlow( std::map< ctr::PDEType, tk::ctr::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the compressible flow system of PDEs
//! \param[inout] cnt std::map of counters for all PDE types
//! \return vector of string pairs describing the PDE configuration
// *****************************************************************************
{
  using eq = newtag::compflow;
  using tk::parameter;
  using tk::parameters;

  auto c = ++cnt[ ctr::PDEType::COMPFLOW ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::COMPFLOW ), "" );

  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< newtag::depvar >()[c] ) );

  nfo.emplace_back( "physics", ctr::Physics().name(
    g_inputdeck.get< newtag::physics >() ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< eq, newtag::problem >() ) );

  auto ncomp = g_inputdeck.get< newtag::ncomp >();
  nfo.emplace_back( "number of components", parameter( ncomp ) );

  const auto scheme = g_inputdeck.get< newtag::scheme >();
  if (scheme != ctr::SchemeType::DiagCG && scheme != ctr::SchemeType::ALECG
    && scheme != ctr::SchemeType::OversetFE)
    nfo.emplace_back( "flux", ctr::Flux().name(
      g_inputdeck.get< newtag::flux >() ) );

  const auto& meshes =
    g_inputdeck.get< newtag::mesh >();

  if (meshes.size() > c) {
    for (const auto& m : meshes) {
      nfo.emplace_back( "mesh", m.get< newtag::filename >() );
    }
  }

  // ICs

  const auto& ic = g_inputdeck.get< newtag::ic >();

  const auto& icbox = ic.get< newtag::box >();
  if (!icbox.empty()) {
    std::size_t bcnt = 0;
    for (const auto& b : icbox) {   // for all boxes configured for this eq
      std::vector< tk::real > box
        { b.get< newtag::xmin >(), b.get< newtag::xmax >(),
          b.get< newtag::ymin >(), b.get< newtag::ymax >(),
          b.get< newtag::zmin >(), b.get< newtag::zmax >() };

      std::string boxname = "IC box " + parameter(bcnt);
      nfo.emplace_back( boxname, parameters( box ) );

      nfo.emplace_back( boxname + " orientation",
        parameters(b.get< newtag::orientation >()) );

      const auto& initiate = b.get< newtag::initiate >();
      auto opt = ctr::Initiate();
      nfo.emplace_back( boxname + ' ' + opt.group(), opt.name(initiate) );

      ++bcnt;
    }
  }

  const auto& icblock = ic.get< newtag::meshblock >();
  for (const auto& b : icblock) {   // for all blocks configured for eq
    std::string blockname = "IC mesh block " +
      parameter(b.get< newtag::blockid >());

    const auto& initiate = b.get< newtag::initiate >();
    auto opt = ctr::Initiate();
    nfo.emplace_back( blockname + ' ' + opt.group(), opt.name(initiate) );
  }

  // BCs

  const auto& bc = g_inputdeck.get< newtag::bc >();
  for (const auto& ib : bc) {
    const auto& stag = ib.get< newtag::stag_point >();
    const auto& radius = ib.get< newtag::radius >();
    if (!stag.empty()) {
      nfo.emplace_back( "Stagnation point(s)", parameters( stag ) );
      nfo.emplace_back( "Stagnation point(s) radii", parameter( radius ) );
    }

    const auto& fs = ib.get< newtag::farfield >();
    if (!fs.empty())
      nfo.emplace_back( "Farfield BC sideset(s)", parameters( fs ) );

    const auto& sym = ib.get< newtag::symmetry >();
    if (!sym.empty())
      nfo.emplace_back( "Symmetry BC sideset(s)", parameters( sym ) );

    const auto& sponge = ib.get< newtag::sponge, newtag::sideset >();
    if (!sponge.empty())
      nfo.emplace_back( "Sponge sideset(s)", parameters( sponge ) );

    const auto& dir = ib.get< newtag::dirichlet >();
    if (!dir.empty())
      nfo.emplace_back( "Dirichlet BC sideset(s)", parameters( dir ) );

    const auto& timedep = ib.get< newtag::timedep >();
    if (!timedep.empty()) {
      for (const auto& bndry : timedep) {
        nfo.emplace_back( "Time dependent BC sideset(s)",
          parameters(bndry.get< newtag::sideset >()) );
      }
    }
  }

  return nfo;
}

void
assignCompFlowGetVars( const std::string& name, tk::GetVarFn& f )
// *****************************************************************************
// Assign functions that compute physics variables from the numerical solution
// for CompFlow
//! \param[in] name Name of variable whose tk::GetVarFn is to be assigned
//! \param[in,out] f Function assigned
// *****************************************************************************
{
  using namespace kw;
  using namespace compflow;

  assign< outvar_density >( name, densityOutVar, f );
  assign< outvar_xvelocity >( name, velocityOutVar<0>, f );
  assign< outvar_yvelocity >( name, velocityOutVar<1>, f );
  assign< outvar_zvelocity >( name, velocityOutVar<2>, f );
  assign< outvar_specific_total_energy >( name, specificTotalEnergyOutVar, f );
  assign< outvar_volumetric_total_energy >
        ( name, volumetricTotalEnergyOutVar, f );
  assign< outvar_xmomentum >( name, momentumOutVar<0>, f );
  assign< outvar_ymomentum >( name, momentumOutVar<1>, f );
  assign< outvar_zmomentum >( name, momentumOutVar<2>, f );
  assign< outvar_pressure >( name, pressureOutVar, f );
}

}  // inciter::
