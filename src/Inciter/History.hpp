// *****************************************************************************
/*!
  \file      src/Inciter/History.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types for collecting history output
  \details   Types for collecting history output.
*/
// *****************************************************************************
#ifndef History_h
#define History_h

#include "Tags.hpp"

namespace inciter {

//! History point data
using HistData = tk::TaggedTuple< brigand::list<
    tag::elem,  std::size_t               //!< Host elem id
  , tag::coord, std::array< tk::real, 3 > //!< Point coordinates
  , tag::fn,    std::array< tk::real, 4 > //!< Shapefunctions evaluated at point
> >;

} // inciter::

#endif // History_h
