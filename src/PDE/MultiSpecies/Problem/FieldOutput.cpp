// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem/FieldOutput.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for multi-species equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible multi-species equations.
*/
// *****************************************************************************
#include "FieldOutput.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "Vector.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "ConfigureMultiSpecies.hpp"
#include "Integrate/Basis.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

std::map< std::string, tk::GetVarFn > MultiSpeciesOutVarFn()
// *****************************************************************************
// Return a map that associates user-specified strings to functions
//! \return Map that associates user-specified strings to functions that compute
//!   relevant quantities to be output to file
// *****************************************************************************
{
  std::map< std::string, tk::GetVarFn > OutFnMap;

  // Allowed strings for user-def field output vars
  OutFnMap["density"] = multispecies::mixDensityOutVar;
  OutFnMap["temperature"] = multispecies::temperatureOutVar;
  OutFnMap["pressure"] = multispecies::pressureOutVar;
  OutFnMap["specific_total_energy"] = multispecies::specificTotalEnergyOutVar;
  OutFnMap["x-velocity"] = multispecies::velocityOutVar<0>;
  OutFnMap["y-velocity"] = multispecies::velocityOutVar<1>;
  OutFnMap["z-velocity"] = multispecies::velocityOutVar<2>;
  for (std::size_t k=0; k<g_inputdeck.get< tag::multispecies, tag::nspec >();
       ++k) {
    const auto mfrac = "mass_fraction_" + std::to_string(k+1);
    OutFnMap[mfrac] = multispecies::massFractionOutVar( k );
  }

  return OutFnMap;
}

std::vector< std::string >
MultiSpeciesFieldNames( std::size_t nspec )
// *****************************************************************************
// Return multi-species field names to be output to file
//! \param[in] nspec Number of species in system
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  for (std::size_t k=0; k<nspec; ++k)
    n.push_back( "density"+std::to_string(k+1)+"_numerical" );
  n.push_back( "x-velocity_numerical" );
  n.push_back( "y-velocity_numerical" );
  n.push_back( "z-velocity_numerical" );
  n.push_back( "pressure_numerical" );
  n.push_back( "soundspeed" );
  n.push_back( "total_energy_density_numerical" );
  n.push_back( "timestep" );

  return n;
}

std::vector< std::string > MultiSpeciesSurfNames()
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
  n.push_back( "wall_normal-temperature_gradient" );

  return n;
}

std::vector< std::vector< tk::real > >
MultiSpeciesSurfOutput(
  const std::size_t nspec,
  const std::size_t rdof,
  const FaceData& fd,
  const tk::Fields& geoFace,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& U,
  const tk::Fields& P )
// *****************************************************************************
//  Return element surface field output (on triangle faces) going to file
//! \param[in] nspec Number of species in this PDE system
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \return Vector of vectors of solution on side set faces to be output to file
// *****************************************************************************
{
  std::vector< std::vector< tk::real > > out;

  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();
  std::vector< EOS > mat_blk;
  initializeSpeciesEoS( mat_blk );
  std::vector< tk::real > ugp( nspec+4, 0.0 );
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];
  std::array< std::vector< tk::real >, 3 > dBdx;
  for (auto& d : dBdx) d.resize( rdof, 0.0 );

  // extract field output along side sets requested
  for (auto s : g_inputdeck.get< tag::field_output, tag::sideset >()) {
    // get face list for side set requested
    auto b = bface.find(static_cast<int>(s));
    if (b == end(bface)) continue;
    const auto& faces = b->second;
    std::vector< tk::real > surfaceSol( faces.size() );
    auto i = out.size();
    out.insert( end(out), 7, surfaceSol );
    std::size_t j = 0;
    for (auto f : faces) {
      Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );
      std::size_t el = static_cast< std::size_t >(esuf[2*f]);

      // access solutions at boundary element
      tk::real rhob(0.0), rhoE(0.0);
      for (std::size_t k=0; k<nspec; ++k) {
        auto rhok = U(el, multispecies::densityDofIdx(nspec,k,rdof,0));
        rhob += rhok;
        ugp[multispecies::densityIdx(nspec,k)] = rhok;
      }
      rhoE = U(el, multispecies::energyDofIdx(nspec,0,rdof,0));
      Mixture mix( nspec, ugp, mat_blk );

      out[i+0][j] = rhob;
      out[i+1][j] = U(el, multispecies::momentumDofIdx(nspec,0,rdof,0))/rhob;
      out[i+2][j] = U(el, multispecies::momentumDofIdx(nspec,1,rdof,0))/rhob;
      out[i+3][j] = U(el, multispecies::momentumDofIdx(nspec,2,rdof,0))/rhob;
      out[i+4][j] = rhoE;
      out[i+5][j] = mix.pressure( rhob,
        P(el, multispecies::temperatureDofIdx(nspec,0,rdof,0)) );
      // Temperature gradients are only computed for second- and higher-order
      // schemes.
      if (rdof >= 4) {
        std::array< std::array< tk::real, 3 >, 4 > coordel {{
          {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
          {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
          {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
          {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }}
        }};
        auto jacInv =
          tk::inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );
        tk::eval_dBdx_p1( rdof, jacInv, dBdx );
        for (std::size_t idir=0; idir<3; ++idir) {
          auto dTdx =
              dBdx[idir][1] *
                P(el, multispecies::temperatureDofIdx(nspec,0,rdof,1))
            + dBdx[idir][2] *
                P(el, multispecies::temperatureDofIdx(nspec,0,rdof,2))
            + dBdx[idir][3] *
                P(el, multispecies::temperatureDofIdx(nspec,0,rdof,3));
          out[i+6][j] += dTdx * geoFace(f,idir+1);
        }
      }
      ++j;
    }
  }

  return out;
}

std::vector< std::string > MultiSpeciesHistNames()
// *****************************************************************************
// Return time history field names to be output to file
//! \note Every time history point will output these fields.
//! \return Vector of strings labelling time history fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;
  auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

  n.push_back( "density" );
  n.push_back( "x-velocity" );
  n.push_back( "y-velocity" );
  n.push_back( "z-velocity" );
  n.push_back( "energy" );
  n.push_back( "pressure" );
  for (std::size_t k=0; k<nspec; ++k)
    n.push_back( "massfrac"+std::to_string(k+1) );

  return n;
}

std::vector< std::string > MultiSpeciesDiagNames(std::size_t nspec)
// *****************************************************************************
// Return diagnostic var names to be output to file
//! \param[in] nspec Number of species in systen
//! \return Vector of strings labelling diagnostic fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  for (std::size_t k=0; k<nspec; ++k)
    n.push_back( "fr"+std::to_string(k+1) );
  n.push_back( "ru" );
  n.push_back( "rv" );
  n.push_back( "rw" );
  n.push_back( "re0" );

  return n;
}

} //inciter::
