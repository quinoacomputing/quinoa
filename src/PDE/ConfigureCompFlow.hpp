// *****************************************************************************
/*!
  \file      src/PDE/ConfigureCompFlow.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
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

namespace inciter {

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
assignCompFlowOutVar( const std::string& name, tk::GetVarFn& f );

/** @name Functions that compute physics variables from the numerical solution for CompFlow */
///@{

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-function"
#endif

//! Compute density for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \return Fluid density ready to be output to file
static tk::GetVarFn::result_type
densityOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset ) {
  return U.extract( 0, offset );
}

//! Compute velocity component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \return Velocity component ready to be output to file
template< tk::ctr::ncomp_t dir >
tk::GetVarFn::result_type
velocityOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset ) {
  using tk::operator/=;
  auto r = U.extract( 0, offset ), u = U.extract( dir, offset );
  u /= r;
  return u;
}

//! Compute total specific energy for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \return Total specific energy ready to be output to file
static tk::GetVarFn::result_type
energyOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset ) {
  return U.extract( 4, offset );
}

//! Compute momentum component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \return Momentum component ready to be output to file
template< tk::ctr::ncomp_t dir >
tk::GetVarFn::result_type
momentumOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset ) {
  return U.extract( dir+1, offset );
}

//! Compute pressure for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the CompFlow
//!   equation system among other systems
//! \return Pressure ready to be output to file
static tk::GetVarFn::result_type
pressureOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset ) {
  auto r = U.extract( 0, offset ), u = U.extract( 1, offset ),
       v = U.extract( 2, offset ), w = U.extract( 3, offset ),
       re = U.extract( 4, offset );
  auto p = r;
  //for (std::size_t i=0; i<U.nunk(); ++i) {
  //  p[i] = eos_pressure<tag::compflow>( system, r[i], u[i], v[i], w[i], re[i] );
  //}
  return p;
}

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//@}

} // inciter::

#endif // ConfigureCompFlow_h
