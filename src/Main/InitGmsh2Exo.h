//******************************************************************************
/*!
  \file      src/Main/InitGmsh2Exo.h
  \author    J. Bakosi
  \date      Wed Mar 19 08:59:23 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh2Exo-specific initialization for main
  \details   Gmsh2Exo-specific initialization for main
*/
//******************************************************************************
#ifndef InitGmsh2Exo_h
#define InitGmsh2Exo_h

#include <Print.h>

namespace gmsh2exo {

//! Echo TPL version information
void echoTPL(const tk::Print& print);

} // quinoa::

#endif // InitGmsh2Exo_h
