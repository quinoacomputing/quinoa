// *****************************************************************************
/*!
  \file      src/Base/Table.h
  \author    J. Bakosi
  \date      Thu 15 Dec 2016 12:41:19 PM MST
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Type alias for declaring, defining, and storing a discrete y = f(x)
             function
  \details   Type alias for declaring, defining, and storing a discrete y = f(x)
             function.
*/
// *****************************************************************************
#ifndef Table_h
#define Table_h

#include <vector>
#include <utility>

#include "Types.h"

namespace tk {

//! Type alias for declaring, defining, and storing a discrete y = f(x) function
using Table = std::vector< std::pair< tk::real, tk::real > >;

} // tk::

#endif // Table_h
