// *****************************************************************************
/*!
  \file      src/Base/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Toolkit-level type definitions
  \details   Toolkit-level type definitions
*/
// *****************************************************************************
#ifndef Types_h
#define Types_h

#include <cstddef>

namespace tk {

//! Real number type used throughout the whole code.
// TODO Test with single precision and possibly others.
using real = double;

//! uint type used throughout the whole code.
using ncomp_t = std::size_t;

} // tk::

#endif // Types_h
