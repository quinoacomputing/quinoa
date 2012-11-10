//******************************************************************************
/*!
  \file      src/Base/QuinoaTypes.h
  \author    J. Bakosi
  \date      Fri 09 Nov 2012 07:26:11 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Code-global type definitions
  \details   Code-global type definitions
*/
//******************************************************************************
#ifndef QuinoaTypes_h
#define QuinoaTypes_h

namespace Quinoa {

using real = double;

// // disregard Intel compiler's remarks
// #ifdef __INTEL_COMPILER
// // operands are evaluated in unspecified order
// #pragma warning(disable:981)
// // value copied to temporary, reference to temporary used
// #pragma warning(disable:383)
// #endif

} // namespace Quinoa

#endif // QuinoaTypes_h
