// *****************************************************************************
/*!
  \file      src/Base/Table.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Basic functionality for storing and sampling a discrete y = f(x)
             function
  \details   Basic functionality for storing and sampling a discrete y = f(x)
             function.
*/
// *****************************************************************************
#ifndef Table_h
#define Table_h

#include <vector>
#include <utility>

#include "Types.hpp"

namespace tk {

//! Type alias for storing a discrete y = f(x) function
using Table = std::vector< std::tuple< tk::real, tk::real > >;

//! Type alias for storing a discrete (y1,y2,y3) = f(x) function
using Table3 =
  std::vector< std::tuple< tk::real, tk::real, tk::real, tk::real > >;

//! Sample a discrete y = f(x) function at x
tk::real sample( tk::real x, const tk::Table& table );

} // tk::

#endif // Table_h
