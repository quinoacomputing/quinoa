// *****************************************************************************
/*!
  \file      src/PDE/DGPDE.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions common to discontinuous Galerkin PDE types
  \details   Functions common to discontinuous Galerkin PDE types.
*/
// *****************************************************************************

#include "DGPDE.hpp"
namespace inciter {

void ConfigInletBC( BCStateFn& s,
                          const tk::StateFn& f,
                          const tk::StateFn& gf )
// *****************************************************************************
// Extract information from input deck on inlet boundary conditions, and append
// to BC configuration vector
//! \param[inout] s BC configuration vector
//! \param[in] f Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] gf Function to evaluate the left and right solution gradients
//!   at boundaries
// *****************************************************************************
{
  const auto& bc = g_inputdeck.get< tag::bc >();
  std::vector< std::size_t > v;
  for (const auto& ib : bc) {
    const auto& in = ib.get< tag::inlet >();
    if (!in.empty()) {
      for (const auto& bndry : in) {
        const auto& sideset = bndry.get< tag::sideset >();
        v.insert(v.end(), sideset.begin(), sideset.end());
        s.push_back( { v, f, gf } );
      }
    }
  }
}

[[noreturn]] tk::StateFn::result_type
invalidBC( ncomp_t, const std::vector< EOS >&,
           const std::vector< tk::real >&, tk::real, tk::real, tk::real,
           tk::real, const std::array< tk::real, 3>& )
// *****************************************************************************
//! State function for invalid/un-configured boundary conditions
//! \note The function signature must follow tk::StateFn
// *****************************************************************************
{
  Throw( "Invalid boundary condition set up in input file or the PDE does not "
          "support this BC type" );
}

}  // inciter::
