//******************************************************************************
/*!
  \file      src/Base/Types.h
  \author    J. Bakosi
  \date      Mon Mar 24 13:57:33 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Toolkit-level type definitions
  \details   Toolkit-level type definitions
*/
//******************************************************************************
#ifndef Types_h
#define Types_h

#include <array>

namespace tk {

using real = double;
using point = std::array< tk::real, 3 >;

} // tk::

#endif // Types_h
