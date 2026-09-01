// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/FieldOutput.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for multi-material equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible multi-material equations.
*/
// *****************************************************************************
#include "FieldOutput.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "Vector.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "ConfigureMultiMat.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

std::map< std::string, tk::GetVarFn > MultiMatOutVarFn()
// *****************************************************************************
// Return a map that associates user-specified strings to functions
//! \return Map that associates user-specified strings to functions that compute
//!   relevant quantities to be output to file
// *****************************************************************************
{
  std::map< std::string, tk::GetVarFn > OutFnMap;

  // Allowed strings for user-def field output vars
  OutFnMap["density"] = multimat::bulkDensityOutVar;
  OutFnMap["pressure"] = multimat::bulkPressureOutVar;
  OutFnMap["specific_total_energy"] = multimat::bulkSpecificTotalEnergyOutVar;
  OutFnMap["x-velocity"] = multimat::velocityOutVar<0>;
  OutFnMap["y-velocity"] = multimat::velocityOutVar<1>;
  OutFnMap["z-velocity"] = multimat::velocityOutVar<2>;
  OutFnMap["material_indicator"] = multimat::matIndicatorOutVar;
  // Cauchy stress tensor
  OutFnMap["stress11"] = multimat::stressOutVar<0,0>;
  OutFnMap["stress12"] = multimat::stressOutVar<0,1>;
  OutFnMap["stress13"] = multimat::stressOutVar<0,2>;
  OutFnMap["stress21"] = multimat::stressOutVar<1,0>;
  OutFnMap["stress22"] = multimat::stressOutVar<1,1>;
  OutFnMap["stress23"] = multimat::stressOutVar<1,2>;
  OutFnMap["stress31"] = multimat::stressOutVar<2,0>;
  OutFnMap["stress32"] = multimat::stressOutVar<2,1>;
  OutFnMap["stress33"] = multimat::stressOutVar<2,2>;
  // Inverse deformation gradient tensor
  OutFnMap["g11"] = multimat::defGradOutVar<0,0>;
  OutFnMap["g12"] = multimat::defGradOutVar<0,1>;
  OutFnMap["g13"] = multimat::defGradOutVar<0,2>;
  OutFnMap["g21"] = multimat::defGradOutVar<1,0>;
  OutFnMap["g22"] = multimat::defGradOutVar<1,1>;
  OutFnMap["g23"] = multimat::defGradOutVar<1,2>;
  OutFnMap["g31"] = multimat::defGradOutVar<2,0>;
  OutFnMap["g32"] = multimat::defGradOutVar<2,1>;
  OutFnMap["g33"] = multimat::defGradOutVar<2,2>;
  // Damage
  OutFnMap["damage"] = multimat::damageOutVar;

  return OutFnMap;
}

std::vector< std::string > MultiMatSurfNames()
// *****************************************************************************
//  Return surface field names to be output to file
//! \note Every surface will output these fields.
//! \return Vector of strings labelling surface fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  n.push_back( "density" );
  n.push_back( "x-velocity" );
  n.push_back( "y-velocity" );
  n.push_back( "z-velocity" );
  n.push_back( "specific_total_energy" );
  n.push_back( "pressure" );

  return n;
}

std::vector< std::vector< tk::real > >
MultiMatSurfOutput(
  const std::size_t nmat,
  const std::size_t rdof,
  const FaceData& fd,
  const tk::Fields& /*geoFace*/,
  const std::vector< std::size_t >& /*inpoel*/,
  const tk::UnsMesh::Coords& /*coord*/,
  const tk::Fields& U,
  const tk::Fields& P )
// *****************************************************************************
//  Return element surface field output (on triangle faces) going to file
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] fd Face connectivity and boundary conditions object
// //! \param[in] geoFace Face geometry array
// //! \param[in] inpoel Element-node connectivity
// //! \param[in] coord Array of nodal coordinates
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \return Vector of vectors of solution on side set faces to be output to file
// *****************************************************************************
{
  std::vector< std::vector< tk::real > > out;

  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();

  // extract field output along side sets requested
  for (auto s : g_inputdeck.get< tag::field_output, tag::sideset >()) {
    // get face list for side set requested
    auto b = bface.find(static_cast<int>(s));
    if (b == end(bface)) continue;
    const auto& faces = b->second;
    std::vector< tk::real > surfaceSol( faces.size() );
    auto i = out.size();
    out.insert( end(out), 6, surfaceSol );
    std::size_t j = 0;
    for (auto f : faces) {
      Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );
      std::size_t el = static_cast< std::size_t >(esuf[2*f]);

      // access solutions at boundary element
      tk::real rhob(0.0), rhoE(0.0), pb(0.0);
      for (std::size_t k=0; k<nmat; ++k) {
        rhob += U(el, densityDofIdx(nmat,k,rdof,0));
        rhoE += U(el, energyDofIdx(nmat,k,rdof,0));
        pb += P(el, pressureDofIdx(nmat,k,rdof,0));
      }

      out[i+0][j] = rhob;
      out[i+1][j] = P(el, velocityDofIdx(nmat,0,rdof,0));
      out[i+2][j] = P(el, velocityDofIdx(nmat,1,rdof,0));
      out[i+3][j] = P(el, velocityDofIdx(nmat,2,rdof,0));
      out[i+4][j] = rhoE;
      out[i+5][j] = pb;
      ++j;
    }
  }

  return out;
}

std::vector< std::string > MultiMatHistNames()
// *****************************************************************************
// Return time history field names to be output to file
//! \note Every time history point will output these fields.
//! \return Vector of strings labelling time history fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

  n.push_back( "density" );
  n.push_back( "x-velocity" );
  n.push_back( "y-velocity" );
  n.push_back( "z-velocity" );
  n.push_back( "energy" );
  n.push_back( "pressure" );
  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "volfrac"+std::to_string(k+1) );

  return n;
}

std::vector< std::string > MultiMatDiagNames(std::size_t nmat)
// *****************************************************************************
// Return diagnostic var names to be output to file
//! \param[in] nmat Number of materials in systen
//! \return Vector of strings labelling diagnostic fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;
  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "f"+std::to_string(k+1) );
  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "fr"+std::to_string(k+1) );
  n.push_back( "fru" );
  n.push_back( "frv" );
  n.push_back( "frw" );
  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "fre"+std::to_string(k+1) );
  for (std::size_t k=0; k<nmat; ++k) {
    if (solidx[k]) {
      for (std::size_t i=1; i<=3; ++i)
        for (std::size_t j=1; j<=3; ++j)
          n.push_back( "g"+std::to_string(k+1)+
            "_"+std::to_string(i)+std::to_string(j) );
    }
  }

  return n;
}

} //inciter::
