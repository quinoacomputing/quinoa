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
infoMultiMat( std::map< ctr::PDEType, tk::ncomp_t >& cnt )
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

  nfo.emplace_back( "physics", ctr::Physics().name(
    g_inputdeck.get< eq, tag::physics >() ) );

  nfo.emplace_back( "problem", ctr::Problem().name(
    g_inputdeck.get< eq, tag::problem >() ) );

  nfo.emplace_back( "flux", ctr::Flux().name(
    g_inputdeck.get< tag::flux >() ) );

  auto nmat = g_inputdeck.get< eq, tag::nmat >();
  nfo.emplace_back( "number of materials", std::to_string( nmat ) );

  auto prelax = g_inputdeck.get< eq, tag::prelax >();
  nfo.emplace_back( "finite pressure relaxation", std::to_string( prelax ) );

  auto intsharp = g_inputdeck.get< eq, tag::intsharp >();
  nfo.emplace_back( "interface sharpening", std::to_string( intsharp ) );

  auto sos_mass_avg = g_inputdeck.get< eq, tag::sos_mass_avg >();
  nfo.emplace_back( "mass average speed of sound",
                    std::to_string( sos_mass_avg ) );

  auto ncomp = g_inputdeck.get< tag::ncomp >();
  nfo.emplace_back( "number of components", std::to_string( ncomp ) );

  // Material eos output
  const auto& matprop = g_inputdeck.get< tag::material >();
  for (const auto& mtype : matprop) {
    const auto& m_id = mtype.get< tag::id >();
    ctr::Material opt;
    nfo.emplace_back( opt.name( mtype.get< tag::eos >() ),
      std::to_string(m_id.size())+" materials" );
    nfo.emplace_back( "material id", parameters( m_id ) );
  }

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

      nfo.emplace_back( boxname + " material id",
                        parameter( b.get< tag::materialid >() ) );

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
