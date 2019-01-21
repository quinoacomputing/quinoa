// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Basis.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing the Dubiner basis functions in DG methods
  \details   This file contains functionality for computing the Dubiner basis
     functions and relating coordinates transformation functions used in 
     discontinuous Galerkin methods for variaous orders of numerical representation.
*/
// *****************************************************************************

#ifndef Basis_h
#define Basis_h

#include "Types.h"
#include "Vector.h"
#include "Fields.h"
#include "FaceData.h"
#include "UnsMesh.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;
using bcconf_t = kw::sideset::info::expect::type;

//! Compute the determination of Jacobian matrix for the tetrahedron element
tk::real
eval_det ( const std::size_t e,
           const std::vector< tk::real >& cx,
           const std::vector< tk::real >& cy,
           const std::vector< tk::real >& cz,
           const std::vector< std::size_t >& inpoel,
           std::array< std::array< tk::real, 3>, 4 >& coordel );

//! Compute the coordinates of quadrature points in physical domain
void
eval_gp ( const std::size_t igp,
          const std::array< std::array< tk::real, 3>, 3 >& coordfa,
          const std::array< std::vector< tk::real >, 2 >& coordgp,
          std::array < tk::real, 3 >& gp );

//! Compute the coordinates of quadrature points in reference domian
void
eval_xi ( const std::array< std::array< tk::real, 3>, 4 >& coordel,
          const tk::real detT,
          const std::array < tk::real, 3 >& gp,
          tk::real& xi,
          tk::real& eta,
          tk::real& zeta );

//! Compute the Dubiner basis functions
void
eval_basis( const tk::real xi,
            const tk::real eta,
            const tk::real zeta,
            std::array< tk::real, 10>& B );

//! Compute the state variables for the tetrahedron element
void
eval_state ( ncomp_t ncomp,
             ncomp_t offset,
             const std::size_t e,
             const Fields& U,
             const Fields& limFunc,
             std::array< tk::real, 10>& B,
             std::vector< tk::real >& state );

} // tk::

#endif // Basis_h
