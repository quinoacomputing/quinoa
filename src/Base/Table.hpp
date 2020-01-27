// *****************************************************************************
/*!
  \file      src/Base/Table.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
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

//! Type alias for declaring, defining, and storing a discrete y = f(x) function
using Table = std::vector< std::pair< tk::real, tk::real > >;

//! Sample a discrete y = f(x) function at x
tk::real sample( tk::real x, const tk::Table& table );

} // tk::

#endif // Table_h
