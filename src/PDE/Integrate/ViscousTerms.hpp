// *****************************************************************************
/*!
  \file      src/PDE/Integrate/ViscousTerms.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing viscous surface terms of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing viscous surface
     terms of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************
#ifndef ViscousTerms_h
#define ViscousTerms_h

#include <array>
#include <vector>

#include "Fields.hpp"
#include "EoS/EOS.hpp"

namespace tk {

//! Viscous terms for multispecies PDEs using the modified-gradient approach
class MultiSpeciesViscousTermsP0P1 {
  public:
    //! Constructor
    MultiSpeciesViscousTermsP0P1(
      std::size_t nspec,
      const std::vector< inciter::EOS >& mat_blk,
      std::size_t rdof,
      const Fields& U,
      const Fields& P );

    //! Return local polynomial order
    std::size_t localDof( std::size_t ) const { return m_rdof; }

    //! Reconstruct conserved variables and primitives at a face point
    std::vector< tk::real >
    stateAt( std::size_t e,
      std::size_t ndof,
      const std::vector< tk::real >& B ) const;

    //! Compute the multispecies viscous flux at an interior face
    std::vector< tk::real >
    interiorFlux(
      const std::array< std::size_t, 2 >& elem,
      const std::array< tk::real, 3 >& fn,
      const std::array< tk::real, 3 >& gp,
      const std::array< std::array< tk::real, 3 >, 2 >& ref_gp,
      const std::array< std::array< tk::real, 3 >, 2 >& centroid,
      const std::array< std::vector< tk::real >, 2 >& B,
      const std::array< std::array< std::vector< tk::real >, 3 >, 2 >& dBdx )
      const;

  private:
    std::size_t m_nspec;
    const std::vector< inciter::EOS >& m_mat_blk;
    std::size_t m_rdof;
    const Fields& m_U;
    const Fields& m_P;
};

} // tk::

#endif // ViscousTerms_h
