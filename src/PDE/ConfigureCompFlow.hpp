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
#include "EoS/GetMatProp.hpp"

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
//! \return Fluid density ready to be output to file
static tk::GetVarFn::result_type
densityOutVar( const tk::Fields& U, std::size_t )
{
  return U.extract_comp( 0 );
}

//! Compute velocity component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Velocity component ready to be output to file
template< tk::ctr::ncomp_t dir >
tk::GetVarFn::result_type
velocityOutVar( const tk::Fields& U, std::size_t rdof )
{
  using tk::operator/=;
  auto r = U.extract_comp( 0 ), u = U.extract_comp( (dir+1)*rdof );
  u /= r;
  return u;
}

//! Compute volumetric total energy (energy per unit volume) for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Volumetric total energy ready to be output to file
static tk::GetVarFn::result_type
volumetricTotalEnergyOutVar( const tk::Fields& U, std::size_t rdof )
{
  return U.extract_comp( 4*rdof );
}

//! Compute specific total energy (energy per unit mass) for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Specific total energy ready to be output to file
static tk::GetVarFn::result_type
specificTotalEnergyOutVar( const tk::Fields& U, std::size_t rdof )
{
  using tk::operator/=;
  auto r = U.extract_comp( 0 ), e = U.extract_comp( 4*rdof );
  e /= r;
  return e;
}

//! Compute momentum component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Momentum component ready to be output to file
template< tk::ctr::ncomp_t dir >
tk::GetVarFn::result_type
momentumOutVar( const tk::Fields& U, std::size_t rdof )
{
  return U.extract_comp( (dir+1)*rdof );
}

//! Compute pressure for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Pressure ready to be output to file
static tk::GetVarFn::result_type
pressureOutVar( const tk::Fields& U, std::size_t rdof )
{
  using tk::operator/=;
  auto r = U.extract_comp( 0 ),
       u = U.extract_comp( 1*rdof ),
       v = U.extract_comp( 2*rdof ),
       w = U.extract_comp( 3*rdof ),
       re = U.extract_comp( 4*rdof );
  u /= r;
  v /= r;
  w /= r;
  auto p = r;
  auto sys = tk::cref_find( g_inputdeck.get< tag::sys >(), 0 );
  for (std::size_t i=0; i<U.nunk(); ++i) {
    // \brief This uses the old eos_pressure call for now, because we didn't 
    // want to change the GetVarFn function signature right now. It's only in
    // the single material CompFlow class, so it shouldn't need multi-material
    // EOSs anyway.
    auto g = gamma< tag::compflow >(sys);
    auto p_c = pstiff< tag::compflow >(sys);
    p[i] = (re[i] - 0.5 * r[i] * (u[i]*u[i] + v[i]*v[i] + w[i]*w[i]) - p_c)
                            * (g-1.0) - p_c;
//    p[i] = m_mat_blk[0]->eos_pressure( sys, r[i], u[i], v[i], w[i], re[i] );
  }
  return p;
}

} // compflow::

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//@}

} // inciter::

#endif // ConfigureCompFlow_h
