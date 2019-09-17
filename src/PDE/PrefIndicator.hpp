// *****************************************************************************
/*!
  \file      src/PDE/Indicator.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
#include "Inciter/Options/PrefIndicator.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Evaluate the adaptive indicator and mark the ndof for each element
void eval_ndof( std::size_t nunk,
                const tk::UnsMesh::Coords& coord,
                const std::vector< std::size_t >& inpoel,
                const inciter::FaceData& fd,
                const tk::Fields& unk,
                inciter::ctr::PrefIndicatorType indicator,
                std::size_t ndof,
                std::size_t ndofmax,
                tk::real tolref,
                std::vector< std::size_t >& ndofel );

} // inciter::

#endif // Indicator_h
