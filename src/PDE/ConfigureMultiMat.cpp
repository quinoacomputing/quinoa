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
#include "MultiMat/Physics/FV.hpp"
#include "MultiMat/DGMultiMat.hpp"
#include "MultiMat/FVMultiMat.hpp"
#include "MultiMat/Problem.hpp"
#include "Inciter/Options/Material.hpp"

namespace inciter {

void
registerMultiMat( DGFactory& df, FVFactory& ff,
  std::set< ctr::PDEType >& fvt, std::set< ctr::PDEType >& dgt )
// *****************************************************************************
// Register multi-material compressible flow PDE into PDE factory
//! \param[in,out] df Discontinuous Galerkin PDE factory to register to
//! \param[in,out] ff Finite volume PDE factory to register to
//! \param[in,out] dgt Counters for equation types registered into DG factory
//! \param[in,out] fvt Counters for equation types registered into FV factory
// *****************************************************************************
{
  // Construct vector of vectors for all possible policies
  using DGMultiMatPolicies =
    tk::cartesian_product< dg::MultiMatPhysics, MultiMatProblems >;
  // Register PDEs for all combinations of policies
  brigand::for_each< DGMultiMatPolicies >(
    registerDG< dg::MultiMat >( df, dgt, ctr::PDEType::MULTIMAT ) );

  // Construct vector of vectors for all possible policies
  using FVMultiMatPolicies =
    tk::cartesian_product< fv::MultiMatPhysics, MultiMatProblems >;
  // Register PDEs for all combinations of policies
  brigand::for_each< FVMultiMatPolicies >(
    registerFV< fv::MultiMat >( ff, fvt, ctr::PDEType::MULTIMAT ) );
}

std::vector< std::pair< std::string, std::string > >
infoMultiMat( std::map< ctr::PDEType, tk::ctr::ncomp_t >& cnt )
// *****************************************************************************
//  Return information on the compressible flow system of PDEs
//! \param[inout] cnt std::map of counters for all PDE types
//! \return vector of string pairs describing the PDE configuration
// *****************************************************************************
{
  using eq = newtag::multimat;
  using tk::parameter;
  using tk::parameters;

  auto c = ++cnt[ ctr::PDEType::MULTIMAT ];       // count eqs
  --c;  // used to index vectors starting with 0

  std::vector< std::pair< std::string, std::string > > nfo;

  nfo.emplace_back( ctr::PDE().name( ctr::PDEType::MULTIMAT ), "" );

  nfo.emplace_back( "dependent variable", std::string( 1,
    g_newinputdeck.get< newtag::depvar >()[c] ) );

  nfo.emplace_back( "physics", ctr::Physics().name(
    g_newinputdeck.get< newtag::physics >() ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_newinputdeck.get< eq, newtag::problem >() ) );

  nfo.emplace_back( "flux", ctr::Flux().name(
    g_newinputdeck.get< newtag::flux >() ) );

  auto nmat = g_newinputdeck.get< eq, newtag::nmat >();
  nfo.emplace_back( "number of materials", std::to_string( nmat ) );

  auto prelax = g_newinputdeck.get< eq, newtag::prelax >();
  nfo.emplace_back( "finite pressure relaxation", std::to_string( prelax ) );

  auto intsharp = g_newinputdeck.get< eq, newtag::intsharp >();
  nfo.emplace_back( "interface sharpening", std::to_string( intsharp ) );

  auto ncomp = g_newinputdeck.get< newtag::ncomp >();
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  // Material eos output
  const auto& matprop = g_newinputdeck.get< newtag::material >();
  for (const auto& mtype : matprop) {
    const auto& m_id = mtype.get< newtag::id >();
    ctr::Material opt;
    nfo.emplace_back( opt.name( mtype.get< newtag::eos >() ),
      std::to_string(m_id.size())+" materials" );
    nfo.emplace_back( "material id", parameters( m_id ) );
  }

  // ICs and IC-boxes

  const auto& ic = g_newinputdeck.get< newtag::ic >();

  const auto& bgmatidic = ic.get< newtag::materialid >();
  nfo.emplace_back( "IC background material id", parameter( bgmatidic ) );

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

      nfo.emplace_back( boxname + " material id",
                        parameter( b.get< newtag::materialid >() ) );

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

    nfo.emplace_back( blockname + " material id",
                      parameter( b.get< newtag::materialid >() ) );
    const auto& initiate = b.get< newtag::initiate >();
    auto opt = ctr::Initiate();
    nfo.emplace_back( blockname + ' ' + opt.group(), opt.name(initiate) );
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
