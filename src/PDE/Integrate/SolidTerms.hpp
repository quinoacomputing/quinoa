// *****************************************************************************
/*!
  \file      src/PDE/Integrate/SolidTerms.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing integrals for the RHS terms
             of the inverse deformation equations in the multi-material
             solid dynamics solver
  \details   This file contains routines to integrate right hand side terms
             to the inverse deformation equations for the multi-material
             solid dynamics solver.
*/
// *****************************************************************************
#ifndef SolidTerms_h
#define SolidTerms_h

#include "Types.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "EoS/EOS.hpp"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute Solid Terms with volume integrals
void
solidTermsVolInt( ncomp_t system,
                  std::size_t nmat,
                  const std::vector< inciter::EOS >& mat_blk,
                  const std::size_t ndof,
                  const std::size_t rdof,
                  const std::size_t nelem,
                  const std::vector< std::size_t >& inpoel,
                  const UnsMesh::Coords& coord,
                  const Fields& geoElem,
                  const Fields& U,
                  const Fields& P,
                  const std::vector< std::size_t >& ndofel,
                  const tk::real dt,
                  Fields& R,
                  int intcompr=0 );

//! Compute Solid Terms with surface integrals
void
solidTermsSurfInt( std::size_t nmat,
                   const std::size_t ndof,
                   const std::size_t rdof,
                   const std::array< tk::real, 3 >& fn,
                   const std::size_t el,
                   const std::size_t er,
                   const std::vector< std::size_t >& solidx,
                   const Fields& geoElem,
                   const Fields& U,
                   const std::array< std::array< tk::real, 3>, 4 > coordel_l,
                   const std::array< std::array< tk::real, 3>, 4 > coordel_r,
                   const std::size_t igp,
                   const std::array< std::vector< tk::real >, 2 >& coordgp,
                   const tk::real dt,
                   std::vector< tk::real >& fl );

// Update the rhs by adding volume integration terms
void
update_rhs( std::size_t nmat,
            const std::size_t ndof,
            const tk::real wt,
            const std::size_t e,
            const std::size_t solidx,
            const std::vector< real >& s,
            Fields& R );

} // tk::

#endif // SolidTerms_h
