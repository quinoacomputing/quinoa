// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for single-material equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible single-material equations.
*/
// *****************************************************************************
#ifndef CompFlowFieldOutput_h
#define CompFlowFieldOutput_h

#include "Fields.hpp"
#include "EoS/EOS.hpp"
#include "History.hpp"
#include "FunctionPrototypes.hpp"
#include "EoS/GetMatProp.hpp"
#include "ContainerUtil.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Return a map that associates user-specified strings to functions
std::map< std::string, tk::GetVarFn > CompFlowOutVarFn();

//! Return surface field names to be output to file
std::vector< std::string > CompFlowSurfNames();

//! Return surface field output going to file
std::vector< std::vector< tk::real > >
CompFlowSurfOutput( const std::vector< EOS >& mat_blk,
                    const std::map< int, std::vector< std::size_t > >& bnd,
                    const tk::Fields& U );

//! Return element surface field output (on triangle faces) going to file
std::vector< std::vector< tk::real > >
CompFlowElemSurfOutput(
  const std::vector< EOS >& mat_blk,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel,
  const tk::Fields& U );

//! Return time history field names to be output to file
std::vector< std::string > CompFlowHistNames();

//! Return time history field output evaluated at time history points
std::vector< std::vector< tk::real > >
CompFlowHistOutput( const std::vector< EOS >& mat_blk,
                    const std::vector< HistData >& h,
                    const std::vector< std::size_t >& inpoel,
                    const tk::Fields& U );

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
template< tk::ncomp_t dir >
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
template< tk::ncomp_t dir >
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
  for (std::size_t i=0; i<U.nunk(); ++i) {
    // \brief This uses the old eos_pressure call for now, because we didn't 
    // want to change the GetVarFn function signature right now. It's only in
    // the single material CompFlow class, so it shouldn't need multi-material
    // EOSs anyway.
    auto g = getmatprop< tag::gamma >();
    auto p_c = getmatprop< tag::pstiff >();
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

} //inciter::

#endif // CompFlowFieldOutput_h
