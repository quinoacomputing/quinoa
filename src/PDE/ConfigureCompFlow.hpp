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

//! Compute density for output to file
//! \note Must follow the signature in tk::GetVarFn
tk::GetVarFn::result_type
densityOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset );

//@}

} // inciter::

#endif // ConfigureCompFlow_h
