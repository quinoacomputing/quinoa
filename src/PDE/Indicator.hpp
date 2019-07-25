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

#include <array>
#include <vector>

#include "Types.hpp"
#include "Fields.hpp"
#include "DerivedData.hpp"
#include "Integrate/Quadrature.hpp"
#include "Integrate/Basis.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

//! Evaluate the adaptive indicator and mark the ndof for each element
void eval_ndof( const std::size_t nunk,
                const std::vector< int >& esuel,
                const tk::Fields& unk,
                std::vector< std::size_t >& ndofel );

} // inciter::

#endif // Indicator_h
