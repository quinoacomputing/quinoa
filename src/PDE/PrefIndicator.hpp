// *****************************************************************************
/*!
  \file      src/Inciter/PrefIndicator.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Adaptive indicators for p-adaptive discontiunous Galerkin methods
  \details   This file contains functions that provide adaptive indicator
    function calculations for marking the number of degree of freedom of each
    element.
*/
// *****************************************************************************
#ifndef Indicator_h
#define Indicator_h

#include <vector>

#include "Types.hpp"
#include "Fields.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

using ncomp_t = tk::ncomp_t;

//! Evaluate the spectral-decay indicator and mark the ndof for each element
void spectral_decay( std::size_t nmat,
                     std::size_t nunk,
                     const std::vector< int >& esuel,
                     const tk::Fields& unk,
                     const tk::Fields& prim,
                     std::size_t ndof,
                     std::size_t ndofmax,
                     tk::real tolref,
                     std::vector< std::size_t >& ndofel );

//! Evaluate the non-conformity indicator and mark the ndof for each element
void non_conformity( std::size_t nunk,
                     std::size_t Nbfac,
                     const std::vector< std::size_t >& inpoel,
                     const tk::UnsMesh::Coords& coord,
                     const std::vector< int >& esuel,
                     const std::vector< int >& esuf,
                     const std::vector< std::size_t >& inpofa,
                     const tk::Fields& unk,
                     std::size_t ndof,
                     std::size_t ndofmax,
                     std::vector< std::size_t >& ndofel );

//! Evaluate the spectral decay indicator for single-material flow
tk::real evalDiscIndicator_CompFlow( std::size_t e,
                                     ncomp_t ncomp,
                                     const std::size_t ndof,
                                     const std::size_t ndofel,
                                     const tk::Fields& unk );

//! Evaluate the spectral decay indicator for multi-material flow
tk::real evalDiscIndicator_MultiMat( std::size_t e,
                                     std::size_t nmat,
                                     ncomp_t ncomp,
                                     ncomp_t nprim,
                                     const std::size_t ndof,
                                     const std::size_t ndofel,
                                     const tk::Fields& unk,
                                     const tk::Fields& prim );
} // inciter::

#endif // Indicator_h
