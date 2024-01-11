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
  using eq = tag::compflow;
  using tk::parameter;
  using tk::parameters;

  auto c = ++cnt[ ctr::PDEType::COMPFLOW ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::COMPFLOW ), "" );

  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, eq, tag::depvar >()[c] ) );

  nfo.emplace_back( "physics", ctr::Physics().name(
    g_inputdeck.get< tag::param, eq, tag::physics >() ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< tag::param, eq, tag::problem >()[c] ) );

  auto ncomp = g_inputdeck.get< tag::component >().get< eq >()[c];
  nfo.emplace_back( "number of components", parameter( ncomp ) );

  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  if (scheme != ctr::SchemeType::DiagCG && scheme != ctr::SchemeType::ALECG
    && scheme != ctr::SchemeType::OversetFE)
    nfo.emplace_back( "flux", ctr::Flux().name(
      g_inputdeck.get< tag::param, eq, tag::flux >().at(c) ) );

  const auto& meshes =
    g_inputdeck.get< tag::param, eq, tag::mesh, tag::filename >();
  if (meshes.size() > c) nfo.emplace_back( "mesh", meshes[c] );

  // Material property output
  //const auto& matprop = g_inputdeck.get< tag::param, eq, tag::material >()[c][0];
  //const auto& m_id = matprop.get< tag::id >();
  //ctr::Material mopt;
  //nfo.emplace_back( mopt.name( matprop.get< tag::eos >() ),
  //  std::to_string(m_id.size()) );

  //nfo.emplace_back( "ratio of specific heats",
  //  parameters(matprop.get< tag::gamma >()) );
  //const auto& cv = matprop.get< tag::cv >();
  //if (!cv.empty())
  //  nfo.emplace_back( "specific heat at constant volume",
  //    parameters(cv) );
  //const auto& pstiff = matprop.get< tag::pstiff >();
  //if (!pstiff.empty())
  //  nfo.emplace_back( "material stiffness",
  //    parameters(pstiff) );

  //// Viscosity is optional: vector may be empty
  //const auto& mu = matprop.get< tag::mu >();
  //if (!mu.empty())
  //  nfo.emplace_back( "dynamic viscosity", parameters( mu ) );

  //// Heat conductivity is optional: vector may be empty
  //const auto& k = matprop.get< tag::k >();
  //if (!k.empty())
  //  nfo.emplace_back( "heat conductivity", parameters( k ) );

  // ICs

  const auto& ic = g_inputdeck.get< tag::param, eq, tag::ic >();

  const auto& bgdensityic = ic.get< tag::density >();
  if (!bgdensityic.empty())
    nfo.emplace_back( "IC background density",
                      parameters( bgdensityic ) );
  const auto& bgvelocityic = ic.get< tag::velocity >();
  if (!bgvelocityic.empty())
    nfo.emplace_back( "IC background velocity",
                      parameters( bgvelocityic ) );
  const auto& bgpressureic = ic.get< tag::pressure >();
  if (!bgpressureic.empty())
    nfo.emplace_back( "IC background pressure",
                      parameters( bgpressureic ) );
  const auto& bgenergyic = ic.get< tag::energy >();
  if (!bgenergyic.empty())
    nfo.emplace_back( "IC background energy",
                      parameters( bgenergyic ) );
  const auto& bgtemperatureic = ic.get< tag::temperature >();
  if (!bgtemperatureic.empty())
    nfo.emplace_back( "IC background temperature",
                      parameters( bgtemperatureic ) );

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

      nfo.emplace_back( boxname + " density",
                        parameter( b.get< tag::density >() ) );
      nfo.emplace_back( boxname + " velocity",
                        parameters( b.get< tag::velocity >() ) );
      nfo.emplace_back( boxname + " pressure",
                        parameter( b.get< tag::pressure >() ) );
      nfo.emplace_back( boxname + " internal energy per unit mass",
                        parameter( b.get< tag::energy >() ) );
      nfo.emplace_back( boxname + " mass",
                        parameter( b.get< tag::mass >() ) );
      nfo.emplace_back( boxname + " internal energy per unit volume",
                        parameter( b.get< tag::energy_content >() ) );
      nfo.emplace_back( boxname + " temperature",
                        parameter( b.get< tag::temperature >() ) );

      const auto& initiate = b.get< tag::initiate >();
      const auto& inittype = initiate.get< tag::init >();
      auto opt = ctr::Initiate();
      nfo.emplace_back( boxname + ' ' + opt.group(), opt.name(inittype) );
      if (inittype == ctr::InitiateType::LINEAR) {
        nfo.emplace_back( boxname + " initiate linear point(s)",
                          parameters( initiate.get< tag::point >() ) );
        nfo.emplace_back( boxname + " initiate linear front width",
                          parameter( initiate.get< tag::front_width >() ) );
        nfo.emplace_back( boxname + " initiate linear velocity",
                          parameter( initiate.get< tag::velocity >() ) );
      }
      ++bcnt;
    }
  }

  const auto& icblock = ic.get< tag::meshblock >();
  for (const auto& b : icblock) {   // for all blocks configured for eq
    std::string blockname = "IC mesh block " +
      parameter(b.get< tag::blockid >());

    nfo.emplace_back( blockname + " volume",
                      parameter( b.get< tag::volume >() ) );
    nfo.emplace_back( blockname + " density",
                      parameter( b.get< tag::density >() ) );
    nfo.emplace_back( blockname + " velocity",
                      parameters( b.get< tag::velocity >() ) );
    nfo.emplace_back( blockname + " pressure",
                      parameter( b.get< tag::pressure >() ) );
    nfo.emplace_back( blockname + " internal energy per unit mass",
                      parameter( b.get< tag::energy >() ) );
    nfo.emplace_back( blockname + " mass",
                      parameter( b.get< tag::mass >() ) );
    nfo.emplace_back( blockname + " internal energy per unit volume",
                      parameter( b.get< tag::energy_content >() ) );
    nfo.emplace_back( blockname + " temperature",
                      parameter( b.get< tag::temperature >() ) );
    const auto& initiate = b.get< tag::initiate >();
    const auto& inittype = initiate.get< tag::init >();
    auto opt = ctr::Initiate();
    nfo.emplace_back( blockname + ' ' + opt.group(), opt.name(inittype) );
  }

  // BCs

  const auto& stag = g_inputdeck.get< tag::param, eq, tag::stag >();
  const auto& spoint = stag.get< tag::point >();
  if (!spoint.empty())
    nfo.emplace_back( "Stagnation point(s)", parameters( spoint ) );
  const auto& sradius = stag.get< tag::radius >();
  if (!sradius.empty())
    nfo.emplace_back( "Stagnation point(s) radii", parameters( sradius ) );

  const auto& skip = g_inputdeck.get< tag::param, eq, tag::skip >();
  const auto& kpoint = skip.get< tag::point >();
  if (!kpoint.empty())
    nfo.emplace_back( "Skip point(s)", parameters( kpoint ) );
  const auto& kradius = skip.get< tag::radius >();
  if (!kradius.empty())
    nfo.emplace_back( "Skip point(s) radii", parameters( kradius ) );

  const auto& fs =
    g_inputdeck.get< tag::param, eq, tag::bc, tag::bcfarfield >();
  if (fs.size() > c) {
    nfo.emplace_back( "Farfield BC sideset(s)", parameters( fs[c] ) );
      nfo.emplace_back( "Farfield BC density", std::to_string(
        g_inputdeck.get< tag::param, eq, tag::farfield_density >() ) );
    const auto& fu =
      g_inputdeck.get< tag::param, eq, tag::farfield_velocity >();
    if (!fu.empty()) {
      nfo.emplace_back( "Farfield BC velocity", parameters(fu) );
    }
    nfo.emplace_back( "Farfield BC pressure", std::to_string(
      g_inputdeck.get< tag::param, eq, tag::farfield_pressure >() ) );
  }

  const auto& sym =
    g_inputdeck.get< tag::param, eq, tag::bc, tag::bcsym >();
  if (sym.size() > c) {
    nfo.emplace_back( "Symmetry BC sideset(s)", parameters( sym[c] ) );

    const auto& sponge = g_inputdeck.get< tag::param, eq, tag::sponge >();
    const auto& ss = sponge.get< tag::sideset >();
    if (!ss.empty()) nfo.emplace_back( "Sponge sideset(s)", parameters( ss ) );
    const auto& spvel = sponge.get< tag::velocity >();
    if (!spvel.empty())
      nfo.emplace_back( "Sponge velocity parameters", parameters( spvel ) );
    const auto& sppre = sponge.get< tag::pressure >();
    if (!sppre.empty())
      nfo.emplace_back( "Sponge pressure parameters", parameters( sppre ) );
  }

  const auto& dir =
    g_inputdeck.get< tag::param, eq, tag::bc, tag::bcdir >();
  if (dir.size() > c)
    nfo.emplace_back( "Dirichlet BC sideset(s)", parameters( dir[c] ) );

  const auto& timedep = g_inputdeck.get< tag::param, eq, tag::bctimedep >();
  if (!timedep.empty()) {
    for (const auto& bndry : timedep) {
      nfo.emplace_back( "Time dependent BC sideset(s)",
        parameters(bndry.get< tag::sideset >()) );
    }
  }

  // FCT

  auto bool_to_string = [](bool B) -> std::string {
    return B ? "true" : "false";
  };

  const auto fct = g_inputdeck.get< tag::discr, tag::fct >();
  if (scheme == ctr::SchemeType::DiagCG && fct) {
    auto sys = g_inputdeck.get< tag::param, eq, tag::sysfct >();
    nfo.emplace_back( "FCT system character", bool_to_string( sys ) );
    if (sys) {
      const auto& sv = g_inputdeck.get< tag::param, eq, tag::sysfctvar >();
      if (!sv.empty()) {
        nfo.emplace_back( "System-FCT variables", parameters( sv ) );
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
