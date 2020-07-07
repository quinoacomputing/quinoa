// *****************************************************************************
/*!
  \file      src/PDE/ConfigureCompFlow.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
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
  using tk::parameters;

  auto c = ++cnt[ ctr::PDEType::COMPFLOW ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::COMPFLOW ), "" );

  nfo.emplace_back( "dependent variable", std::string( 1,
    g_inputdeck.get< tag::param, eq, tag::depvar >()[c] ) );

  nfo.emplace_back( "physics", ctr::Physics().name(
    g_inputdeck.get< tag::param, eq, tag::physics >()[c] ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< tag::param, eq, tag::problem >()[c] ) );

  auto ncomp = g_inputdeck.get< tag::component >().get< eq >()[c];
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  if (scheme != ctr::SchemeType::DiagCG && scheme != ctr::SchemeType::ALECG)
    nfo.emplace_back( "flux", ctr::Flux().name(
      g_inputdeck.get< tag::param, eq, tag::flux >().at(c) ) );

  nfo.emplace_back( "start offset in unknowns array", std::to_string(
    g_inputdeck.get< tag::component >().offset< eq >(c) ) );

  const auto& gamma = g_inputdeck.get< tag::param, eq, tag::gamma >()[c];
  if (!gamma.empty())
    nfo.emplace_back( "ratio of specific heats", std::to_string( gamma[0] )) ;

  const auto& pstiff = g_inputdeck.get< tag::param, eq, tag::pstiff >()[c];
  if (!pstiff.empty())
    nfo.emplace_back( "material stiffness", std::to_string( pstiff[0] ) );

  // Viscosity is optional: the outer vector may be empty
  const auto& mu = g_inputdeck.get< tag::param, eq, tag::mu >();
  if (mu.size() > c)
    nfo.emplace_back( "dynamic viscosity", std::to_string( mu[c][0] ) );

  const auto& cv = g_inputdeck.get< tag::param, eq, tag::cv >()[c];
  if (!cv.empty())
    nfo.emplace_back( "specific heat at const. volume", std::to_string(cv[0]) );

  // Heat conductivity is optional: the outer vector may be empty
  const auto& k = g_inputdeck.get< tag::param, eq, tag::k >();
  if (k.size() > c)
    nfo.emplace_back( "heat conductivity", std::to_string( k[c][0] ) );

  const auto& npar = g_inputdeck.get< tag::param, eq, tag::npar >();
  if (!npar.empty())
    nfo.emplace_back( "number of tracker particles", parameters( npar ) );

  const auto& alpha = g_inputdeck.get< tag::param, eq, tag::alpha >();
  if (!alpha.empty()) nfo.emplace_back( "coeff alpha", parameters( alpha ) );

  const auto& beta =
    g_inputdeck.get< tag::param, eq, tag::beta >();
  if (!beta.empty())
    nfo.emplace_back( "coeff beta", parameters( beta ) );

  const auto& bx = g_inputdeck.get< tag::param, eq, tag::betax >();
  if (!bx.empty()) nfo.emplace_back( "coeff betax", parameters( bx ) );

  const auto& by = g_inputdeck.get< tag::param, eq, tag::betay >();
  if (!by.empty()) nfo.emplace_back( "coeff betay", parameters( by ) );

  const auto& bz = g_inputdeck.get< tag::param, eq, tag::betaz >();
  if (!bz.empty()) nfo.emplace_back( "coeff betaz", parameters( bz ) );

  const auto& r0 = g_inputdeck.get< tag::param, eq, tag::r0 >();
  if (!r0.empty()) nfo.emplace_back( "coeff r0", parameters( r0 ) );

  const auto& ce = g_inputdeck.get< tag::param, eq, tag::ce >();
  if (!ce.empty()) nfo.emplace_back( "coeff ce", parameters( ce ) );

  const auto& kappa = g_inputdeck.get< tag::param, eq, tag::kappa >();
  if (!kappa.empty()) nfo.emplace_back( "coeff k", parameters( kappa ) );

  const auto& p0 = g_inputdeck.get< tag::param, eq, tag::p0 >();
  if (!p0.empty()) nfo.emplace_back( "coeff p0", parameters( p0 ) );

  // ICs

  const auto& ic = g_inputdeck.get< tag::param, eq, tag::ic >();

  const auto& bgdensityic = ic.get< tag::density >();
  if (bgdensityic.size() > c && !bgdensityic[c].empty())
    nfo.emplace_back( "IC background density",
                      std::to_string( bgdensityic[c][0] ) );
  const auto& bgvelocityic = ic.get< tag::velocity >();
  if (bgvelocityic.size() > c && !bgvelocityic[c].empty())
    nfo.emplace_back( "IC background velocity",
                      parameters( bgvelocityic[c] ) );
  const auto& bgpressureic = ic.get< tag::pressure >();
  if (bgpressureic.size() > c && !bgpressureic[c].empty())
    nfo.emplace_back( "IC background pressure",
                      std::to_string( bgpressureic[c][0] ) );

  const auto& icbox = ic.get< tag::box >();
  std::vector< tk::real > box{ icbox.get< tag::xmin >(),
                               icbox.get< tag::xmax >(),
                               icbox.get< tag::ymin >(),
                               icbox.get< tag::ymax >(),
                               icbox.get< tag::zmin >(),
                               icbox.get< tag::zmax >() };
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  if (std::any_of( begin(box), end(box),
        [=]( tk::real p ){ return std::abs(p) > eps; })) {
    nfo.emplace_back( "IC box", parameters( box ) );
  }

  const auto& boxdensityic = icbox.get< tag::density >();
  if (boxdensityic.size() > c)
    nfo.emplace_back( "IC box density",
                      std::to_string( boxdensityic[c][0] ) );
  const auto& boxvelocityic = icbox.get< tag::velocity >();
  if (boxvelocityic.size() > c)
    nfo.emplace_back( "IC box velocity",
                      parameters( boxvelocityic[c] ) );
  const auto& boxpressureic = icbox.get< tag::pressure >();
  if (boxpressureic.size() > c)
    nfo.emplace_back( "IC box pressure",
                      std::to_string( boxpressureic[c][0] ) );
  const auto& boxenergyic = icbox.get< tag::energy >();
  if (boxenergyic.size() > c)
    nfo.emplace_back( "IC box internal energy per unit mass",
                      std::to_string( boxenergyic[c][0] ) );
  const auto& boxmassic = icbox.get< tag::mass >();
  if (boxmassic.size() > c)
    nfo.emplace_back( "IC box mass", std::to_string( boxmassic[c][0] ) );
  const auto& boxenergy_content_ic = icbox.get< tag::energy_content >();
  if (boxenergy_content_ic.size() > c)
    nfo.emplace_back( "IC box internal energy per unit volume",
                      std::to_string( boxenergy_content_ic[c][0] ) );
  const auto& boxtemperatureic = icbox.get< tag::temperature >();
  if (boxtemperatureic.size() > c)
    nfo.emplace_back( "IC box temperature",
                      std::to_string( boxtemperatureic[c][0] ) );

  // BCs

  const auto& bcstag = g_inputdeck.get< tag::param, eq, tag::bcstag >();
  const auto& point = bcstag.get< tag::point >();
  if (point.size() > c)
    nfo.emplace_back( "Stagnation BC point(s)", parameters( point[c] ) );
  const auto& radius = bcstag.get< tag::radius >();
  if (radius.size() > c)
    nfo.emplace_back( "Stagnation BC radii", parameters( radius[c] ) );

  const auto& fs =
    g_inputdeck.get< tag::param, eq, tag::bc, tag::bcfarfield >();
  if (fs.size() > c) {
    nfo.emplace_back( "Farfield BC sideset(s)", parameters( fs[c] ) );
    const auto& fr =
      g_inputdeck.get< tag::param, eq, tag::farfield_density >();
    if (fr.size() > c)
      nfo.emplace_back( "Farfield BC density", std::to_string(fr[c]) );
    const auto& fu =
      g_inputdeck.get< tag::param, eq, tag::farfield_velocity >();
    if (fu.size() > c)
      nfo.emplace_back( "Farfield BC velocity", parameters( fu[c] ) );
    const auto& fp =
      g_inputdeck.get< tag::param, eq, tag::farfield_pressure >();
    if (fp.size() > c)
      nfo.emplace_back( "Farfield BC pressure", std::to_string(fp[c]) );
  }

  const auto& sym =
    g_inputdeck.get< tag::param, eq, tag::bc, tag::bcsym >();
  if (sym.size() > c)
    nfo.emplace_back( "Symmetry BC sideset(s)", parameters( sym[c] ) );

  // FCT

  auto bool_to_string = [](bool b) -> std::string {
    return b ? "true" : "false";
  };

  const auto fct = g_inputdeck.get< tag::discr, tag::fct >();
  if (scheme == ctr::SchemeType::DiagCG && fct) {

    const auto& sys = g_inputdeck.get< tag::param, eq, tag::sysfct >();
    if (sys.size() > c) {
      nfo.emplace_back( "FCT system character", bool_to_string( sys[c] ) );

      if (sys[c]) {     // if system FCT is enabled for this system
        const auto& sv = g_inputdeck.get< tag::param, eq, tag::sysfctvar >();
        if (sv.size() > c) {
          nfo.emplace_back( "System-FCT variables", parameters( sv[c] ) );
        }
      }
    }

  }

  return nfo;
}

}  // inciter::
