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

[[noreturn]] tk::StateFn::result_type
inciter::invalidBC( ncomp_t, ncomp_t, const std::vector< tk::real >&,
           tk::real, tk::real, tk::real, tk::real,
           const std::array< tk::real, 3> &,
           const std::vector< EoS_Base* >& )
// *****************************************************************************
//! State function for invalid/un-configured boundary conditions
//! \note The function signature must follow tk::StateFn
// *****************************************************************************
{
  Throw( "Invalid boundary condition set up in input file or the PDE does not "
          "support this BC type" );
}
