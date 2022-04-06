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
#include "EoS/EoS_Base.hpp"

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
densityOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset, std::size_t )
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
velocityOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset, std::size_t rdof )
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
volumetricTotalEnergyOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset,
                             std::size_t rdof )
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
specificTotalEnergyOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset,
                           std::size_t rdof )
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
momentumOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset, std::size_t rdof )
{
  return U.extract( (dir+1)*rdof, offset );
}

//! Compute pressure for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Pressure ready to be output to file
static tk::GetVarFn::result_type
pressureOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset, std::size_t rdof )
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
    // This uses the old eos_pressure call for now, because we didn't want to 
    // change the GetVarFn function signature right now. It's only in the single
    // material CompFlow class, so it shouldn't need multi-material EOSs anyway.
    p[i] = eos_pressure<tag::compflow>( sys, r[i], u[i], v[i], w[i], re[i] );
//    p[i] = m_mat_blk[0]->eos_pressure( sys, r[i], u[i], v[i], w[i], re[i] );
  return p;
}

} // compflow::

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//@}

} // inciter::

#endif // ConfigureCompFlow_h
