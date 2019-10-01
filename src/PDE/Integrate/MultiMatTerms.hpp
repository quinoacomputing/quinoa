// *****************************************************************************
/*!
  \file      src/PDE/Integrate/MultiMatTerms.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing volume integrals of multi-material terms
     using DG methods
  \details   This file contains functionality for computing volume integrals of
     non-conservative and pressure relaxation terms that appear in the
     multi-material hydrodynamic equations, using the discontinuous Galerkin
     method for various orders of numerical representation.
*/
// *****************************************************************************
#ifndef MultiMatTerms_h
#define MultiMatTerms_h

#include "Basis.hpp"
#include "Types.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute volume integrals of non-conservative terms for multi-material DG
void
nonConservativeInt( ncomp_t system,
                    std::size_t nmat,
                    ncomp_t offset,
                    const std::size_t ndof,
                    const std::size_t rdof,
                    const std::vector< std::size_t >& inpoel,
                    const UnsMesh::Coords& coord,
                    const Fields& geoElem,
                    const Fields& U,
                    const Fields& P,
                    const std::vector< std::vector< tk::real > >& riemannDeriv,
                    const std::vector< std::size_t >& ndofel,
                    Fields& R );

//! Update the rhs by adding the non-conservative term integrals
void
update_rhs_ncn( ncomp_t ncomp,
                ncomp_t offset,
                const std::size_t ndof,
                const std::size_t ndof_el,
                const tk::real wt,
                const std::size_t e,
                const std::array< std::vector<tk::real>, 3 >& dBdx,
                const std::vector< tk::real >& ncf,
                Fields& R );

//! Compute volume integrals of pressure relaxation terms in multi-material DG
void
pressureRelaxationInt( ncomp_t system,
                       ncomp_t ncomp,
                       std::size_t nmat,
                       ncomp_t offset,
                       const std::size_t ndof,
                       const std::size_t rdof,
                       const Fields& geoElem,
                       const Fields& U,
                       const std::vector< std::size_t >& ndofel,
                       const tk::real ct,
                       Fields& R );

} // tk::

#endif // MultiMatTerms_h
