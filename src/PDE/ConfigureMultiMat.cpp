// *****************************************************************************
/*!
  \file      src/PDE/ConfigureMultiMat.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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
#include "ContainerUtil.hpp"
#include "ConfigureMultiMat.hpp"
#include "MultiMat/Physics/DG.hpp"
#include "MultiMat/DGMultiMat.hpp"
#include "MultiMat/Problem.hpp"
#include "InfoMesh.hpp"
#include "Inciter/Options/Material.hpp"

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
infoMultiMat( std::map< ctr::PDEType, tk::ctr::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the compressible flow system of PDEs
//! \param[inout] cnt std::map of counters for all PDE types
//! \return vector of string pairs describing the PDE configuration
// *****************************************************************************
{
  using eq = tag::multimat;
  using tk::parameter;
  using tk::parameters;

  auto c = ++cnt[ ctr::PDEType::MULTIMAT ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::MULTIMAT ), "" );

  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, eq, tag::depvar >()[c] ) );

  infoMesh< eq >( c, nfo );

  nfo.emplace_back( "physics", ctr::Physics().name(
    g_inputdeck.get< tag::param, eq, tag::physics >()[c] ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< tag::param, eq, tag::problem >()[c] ) );

  nfo.emplace_back( "flux", ctr::Flux().name(
    g_inputdeck.get< tag::param, eq, tag::flux >().at(c) ) );

  auto nmat = g_inputdeck.get< tag::param, eq, tag::nmat >()[c];
  nfo.emplace_back( "number of materials", std::to_string( nmat ) );

  auto prelax = g_inputdeck.get< tag::param, eq, tag::prelax >()[c];
  nfo.emplace_back( "finite pressure relaxation", std::to_string( prelax ) );

  if (prelax)
  {
    auto prelax_ts =
      g_inputdeck.get< tag::param, eq, tag::prelax_timescale >()[c];
    nfo.emplace_back( "pressure relaxation time-scale",
                      std::to_string( prelax_ts ) );
  }

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

  nfo.emplace_back( "start offset in unknowns array", std::to_string(
    g_inputdeck.get< tag::component >().offset< eq >(c) ) );

  // Material property output
  const auto& matprop = g_inputdeck.get< tag::param, eq, tag::material >()[c];
  for (const auto& mtype : matprop) {
    const auto& m_id = mtype.get< tag::id >();
    if (mtype.get< tag::eos >() == ctr::MaterialType::STIFFENEDGAS) {
      nfo.emplace_back( "Stiffened Gas EoS materials: ",
        std::to_string(m_id.size()) );
    }
    else if (mtype.get< tag::eos >() == ctr::MaterialType::JWL) {
      nfo.emplace_back( "Jones-Wilkins-Lee (JWL) EoS:",
        std::to_string(m_id.size()) );
    }

    nfo.emplace_back( "material id", parameters( m_id ) );
    nfo.emplace_back( "ratio of specific heats",
      parameters(mtype.get< tag::gamma >()) );
    nfo.emplace_back( "specific heat at constant volume",
      parameters(mtype.get< tag::cv >()) );
    nfo.emplace_back( "material stiffness",
      parameters(mtype.get< tag::pstiff >()) );

    // Viscosity is optional: vector may be empty
    const auto& mu = mtype.get< tag::mu >();
    if (!mu.empty())
      nfo.emplace_back( "dynamic viscosity", parameters( mu ) );

    // Heat conductivity is optional: vector may be empty
    const auto& k = mtype.get< tag::k >();
    if (!k.empty())
      nfo.emplace_back( "heat conductivity", parameters( k ) );
  }

  // ICs and IC-boxes

  const auto& ic = g_inputdeck.get< tag::param, eq, tag::ic >();

  const auto& bgmatidic = ic.get< tag::materialid >();
  if (bgmatidic.size() > c && !bgmatidic[c].empty())
    nfo.emplace_back( "IC background material id",
                      parameter( bgmatidic[c][0] ) );
  const auto& bgdensityic = ic.get< tag::density >();
  if (bgdensityic.size() > c && !bgdensityic[c].empty())
    nfo.emplace_back( "IC background density",
                      parameter( bgdensityic[c][0] ) );
  const auto& bgvelocityic = ic.get< tag::velocity >();
  if (bgvelocityic.size() > c && !bgvelocityic[c].empty())
    nfo.emplace_back( "IC background velocity",
                      parameters( bgvelocityic[c] ) );
  const auto& bgpressureic = ic.get< tag::pressure >();
  if (bgpressureic.size() > c && !bgpressureic[c].empty())
    nfo.emplace_back( "IC background pressure",
                      parameter( bgpressureic[c][0] ) );
  const auto& bgenergyic = ic.get< tag::energy >();
  if (bgenergyic.size() > c && !bgenergyic[c].empty())
    nfo.emplace_back( "IC background energy",
                      parameter( bgenergyic[c][0] ) );
  const auto& bgtemperatureic = ic.get< tag::temperature >();
  if (bgtemperatureic.size() > c && !bgtemperatureic[c].empty())
    nfo.emplace_back( "IC background temperature",
                      parameter( bgtemperatureic[c][0] ) );

  const auto& icbox = ic.get< tag::box >();
  if (icbox.size() > c) {
    std::size_t bcnt = 0;
    for (const auto& b : icbox[c]) {   // for all boxes configured for this eq
      std::vector< tk::real > box
        { b.get< tag::xmin >(), b.get< tag::xmax >(),
          b.get< tag::ymin >(), b.get< tag::ymax >(),
          b.get< tag::zmin >(), b.get< tag::zmax >() };

      std::string boxname = "IC box " + parameter(bcnt);
      nfo.emplace_back( boxname, parameters( box ) );

      nfo.emplace_back( boxname + " material id",
                        parameter( b.get< tag::materialid >() ) );
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

      ++bcnt;
    }
  }

  return nfo;
}

void
assignMultiMatGetVars( const std::string& name, tk::GetVarFn& f )
// *****************************************************************************
// Assign functions that compute physics variables from the numerical solution
// for MultiMat
//! \param[in] name Name of variable whose tk::GetVarFn is to be assigned
//! \param[in,out] f Function assigned
// *****************************************************************************
{
  using namespace kw;
  using namespace multimat;

  assign< outvar_density >( name, bulkDensityOutVar, f );
  assign< outvar_pressure >( name, bulkPressureOutVar, f );
  assign< outvar_specific_total_energy >
        ( name, bulkSpecificTotalEnergyOutVar, f );
  assign< outvar_xvelocity >( name, velocityOutVar<0>, f );
  assign< outvar_yvelocity >( name, velocityOutVar<1>, f );
  assign< outvar_zvelocity >( name, velocityOutVar<2>, f );
  assign< outvar_material_indicator >( name, matIndicatorOutVar, f );
}

}  // inciter::
