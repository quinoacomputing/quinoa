// *****************************************************************************
/*!
  \file      src/PDE/ConfigureCompFlow.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration for compressible flow PDE
  \details   Register and compile configuration for compressible flow PDE.
*/
// *****************************************************************************
#ifndef ConfigureCompFlow_h
#define ConfigureCompFlow_h

#include <set>
#include <map>
#include <vector>

#include "PDEFactory.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/PDE.hpp"
#include "FunctionPrototypes.hpp"
#include "ContainerUtil.hpp"
#include "EoS/EoS.hpp"
#include "Mesh/DerivedData.hpp"
#include "Vector.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Register compressible flow PDEs into PDE factory
void
registerCompFlow( CGFactory& cf,
                  DGFactory& df,
                  std::set< ctr::PDEType >& cgt,
                  std::set< ctr::PDEType >& dgt );

//! Return information on the compressible flow PDE
std::vector< std::pair< std::string, std::string > >
infoCompFlow( std::map< ctr::PDEType, tk::ctr::ncomp_t >& cnt );

//! \brief Assign function that computes physics variables from the
//!   numerical solution for CompFlow
void
assignCompFlowGetVars( const std::string& name, tk::GetVarFn& f );

/** @name Functions that compute physics variables from the numerical solution for CompFlow */
///@{

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-function"
#endif

namespace compflow {

//! Compute density for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \return Fluid density ready to be output to file
static tk::GetVarFn::result_type
densityOutVar( const tk::Fields& U,
               tk::ctr::ncomp_t offset,
               std::size_t,
               const tk::UnsMesh::Coords&,
               const std::vector< std::size_t >&,
               const std::vector< tk::real >& )
{
  return U.extract( 0, offset );
}

//! Compute velocity component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Velocity component ready to be output to file
template< tk::ctr::ncomp_t dir >
tk::GetVarFn::result_type
velocityOutVar( const tk::Fields& U,
                tk::ctr::ncomp_t offset,
                std::size_t rdof,
                const tk::UnsMesh::Coords&,
                const std::vector< std::size_t >&,
                const std::vector< tk::real >& )
{
  using tk::operator/=;
  auto r = U.extract( 0, offset ), u = U.extract( (dir+1)*rdof, offset );
  u /= r;
  return u;
}

//! Compute volumetric total energy (energy per unit volume) for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Volumetric total energy ready to be output to file
static tk::GetVarFn::result_type
volumetricTotalEnergyOutVar( const tk::Fields& U,
                             tk::ctr::ncomp_t offset,
                             std::size_t rdof,
                             const tk::UnsMesh::Coords&,
                             const std::vector< std::size_t >&,
                             const std::vector< tk::real >& )
{
  return U.extract( 4*rdof, offset );
}

//! Compute specific total energy (energy per unit mass) for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Specific total energy ready to be output to file
static tk::GetVarFn::result_type
specificTotalEnergyOutVar( const tk::Fields& U,
                           tk::ctr::ncomp_t offset,
                           std::size_t rdof,
                           const tk::UnsMesh::Coords&,
                           const std::vector< std::size_t >&,
                           const std::vector< tk::real >& )
{
  using tk::operator/=;
  auto r = U.extract( 0, offset ), e = U.extract( 4*rdof, offset );
  e /= r;
  return e;
}

//! Compute momentum component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Momentum component ready to be output to file
template< tk::ctr::ncomp_t dir >
tk::GetVarFn::result_type
momentumOutVar( const tk::Fields& U,
                tk::ctr::ncomp_t offset,
                std::size_t rdof,
                const tk::UnsMesh::Coords&,
                const std::vector< std::size_t >&,
                const std::vector< tk::real >& )
{
  return U.extract( (dir+1)*rdof, offset );
}

//! Compute vorticity component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] vol Nodal volumes
//! \return Momentum component ready to be output to file
template< tk::ctr::ncomp_t dir >
tk::GetVarFn::result_type
vorticityOutVar( const tk::Fields& U,
                 tk::ctr::ncomp_t offset,
                 std::size_t rdof,
                 const tk::UnsMesh::Coords& coord,
                 const std::vector< std::size_t >& inpoel,
                 const std::vector< tk::real >& vol )
{
  using tk::operator/=;
  tk::UnsMesh::Coords vel{ U.extract( 1*rdof, offset ),
                           U.extract( 2*rdof, offset ),
                           U.extract( 3*rdof, offset ) };
  auto r = U.extract( 0, offset );
  vel[0] /= r;
  vel[1] /= r;
  vel[2] /= r;
  auto vort = tk::curl( coord, inpoel, vol, vel );
  return vort[dir];
}

//! Compute length of the vorticity vector for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] vol Nodal volumes
//! \return Momentum component ready to be output to file
static tk::GetVarFn::result_type
vorticityMagOutVar( const tk::Fields& U,
                    tk::ctr::ncomp_t offset,
                    std::size_t rdof,
                    const tk::UnsMesh::Coords& coord,
                    const std::vector< std::size_t >& inpoel,
                    const std::vector< tk::real >& vol )
{
  using tk::operator/=;
  tk::UnsMesh::Coords vel{ U.extract( 1*rdof, offset ),
                           U.extract( 2*rdof, offset ),
                           U.extract( 3*rdof, offset ) };
  auto r = U.extract( 0, offset );
  vel[0] /= r;
  vel[1] /= r;
  vel[2] /= r;
  auto vort = tk::curl( coord, inpoel, vol, vel );
  for (std::size_t i=0; i<U.nunk(); ++i)
    vort[0][i] = tk::length( vort[0][i], vort[1][i], vort[2][i] );
  return vort[0];
}

//! Compute pressure for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Pressure ready to be output to file
static tk::GetVarFn::result_type
pressureOutVar( const tk::Fields& U,
                tk::ctr::ncomp_t offset,
                std::size_t rdof,
                const tk::UnsMesh::Coords&,
                const std::vector< std::size_t >&,
                const std::vector< tk::real >& )
{
  using tk::operator/=;
  auto r = U.extract( 0, offset ),
       u = U.extract( 1*rdof, offset ),
       v = U.extract( 2*rdof, offset ),
       w = U.extract( 3*rdof, offset ),
       re = U.extract( 4*rdof, offset );
  u /= r;
  v /= r;
  w /= r;
  auto p = r;
  auto sys = tk::cref_find( g_inputdeck.get< tag::sys >(), offset );
  for (std::size_t i=0; i<U.nunk(); ++i)
    p[i] = eos_pressure<tag::compflow>( sys, r[i], u[i], v[i], w[i], re[i] );
  return p;
}

} // compflow::

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//@}

} // inciter::

#endif // ConfigureCompFlow_h
