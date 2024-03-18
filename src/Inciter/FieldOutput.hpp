// *****************************************************************************
/*!
  \file      src/Inciter/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Extract field output for inciter
  \details   Extract field output for inciter.
*/
// *****************************************************************************
#ifndef FieldOutput_h
#define FieldOutput_h

#include "Types.hpp"
#include "Fields.hpp"
#include "Centering.hpp"
#include "ContainerUtil.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "UnsMesh.hpp"
#include "Discretization.hpp"

namespace inciter {

extern ctr::New2InputDeck g_newinputdeck;

//! Collect field output names from numerical solution based on user input
std::vector< std::string >
numericFieldNames( tk::Centering c, char depvar = 0 );

//! Collect field output from numerical solution based on user input
std::vector< std::vector< tk::real > >
numericFieldOutput( const tk::Fields& U, tk::Centering c,
                    const tk::Fields& P = tk::Fields(),
                    char depvar = 0 );

//! Evaluate solution on incoming (a potentially refined) mesh
void
evalSolution(
  const Discretization& D,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, std::size_t >& addedTets,
  const std::vector< std::size_t >& ndofel,
  const tk::Fields& U,
  const tk::Fields& P,
  tk::Fields& uElemfields,
  tk::Fields& pElemfields,
  tk::Fields& uNodefields,
  tk::Fields& pNodefields );

//! Collect field output names from analytic solutions based on user input
//! \tparam PDE Partial differential equation type
//! \param[in] eq PDE whose analytic solution field names to query
//! \param[in] c Extract variables only with this centering
//! \param[in,out] f Output field names augmented
template< class PDE >
void
analyticFieldNames( const PDE& eq,
                    tk::Centering c,
                    std::vector< std::string >& f )
{
  for (const auto& v : g_newinputdeck.get< newtag::field_output, newtag::outvar >())
    if (v.centering == c && v.analytic())
      tk::concat( eq.analyticFieldNames(), f );
}

//! Collect field output from analytic solutions based on user input
//! \tparam PDE Partial differential equation type
//! \param[in] eq PDE whose analytic solution to output
//! \param[in] c Extract variables only with this centering
//! \param[in] x x coordinates at which to evaluate the analytic solution
//! \param[in] y y coordinates at which to evaluate the analytic solution
//! \param[in] z z coordinates at which to evaluate the analytic solution
//! \param[in] t Physical time at which to evaluate the analytic solution
//! \param[in,out] f Output fields augmented by analytic solutions requested
template< class PDE >
void
analyticFieldOutput( const PDE& eq,
                     tk::Centering c,
                     const std::vector< tk::real >& x,
                     const std::vector< tk::real >& y,
                     const std::vector< tk::real >& z,
                     tk::real t,
                     std::vector< std::vector< tk::real > >& f )
{
  for (const auto& v : g_newinputdeck.get< newtag::field_output, newtag::outvar >()) {
    if (v.centering == c && v.analytic()) {
      auto ncomp = eq.analyticSolution( x[0], y[0], z[0], t ).size();
      f.resize( f.size() + ncomp, std::vector< tk::real >( x.size() ) );
      for (std::size_t i=0; i<x.size(); ++i) {
        auto s = eq.analyticSolution( x[i], y[i], z[i], t );
        for (std::size_t j=0; j<ncomp; ++j) f[f.size()-ncomp+j][i] = s[j];
      }
    }
  }
}

} // inciter::

#endif // FieldOutput_h
