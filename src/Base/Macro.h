//******************************************************************************
/*!
  \file      src/Base/Macro.h
  \author    J. Bakosi
  \date      Sat 02 Mar 2013 10:37:03 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Macro definitions
  \details   Macro definitions
*/
//******************************************************************************
#ifndef Macro_h
#define Macro_h

namespace Quinoa {

#define SWAP(a,b,tmp) {tmp=a; a=b; b=tmp;}
#define IGNORE(expr) (static_cast<void>(expr))

} // namespace Quinoa

#endif // Macro_h
