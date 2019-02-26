// *****************************************************************************
/*!
  \file      src/UnitTest/QuietCerr.C
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Charm++ nodegroup to quiet std::cerr in a thread-safe fashion
  \details   Charm++ nodegroup to quiet std::cerr in a thread-safe fashion.
*/
// *****************************************************************************

#include <iostream>
#include <sstream>

#include "QuietCerr.h"

namespace tk {

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

//! std::tringstream used to quiet std::cerr's stream by redirecting to it
std::stringstream cerr_quiet;
//! std::streambuf used to store state of std::cerr before redirecting it
std::streambuf* cerr_old;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

}

using tk::QuietCerr;

void
QuietCerr::quiet()
// *****************************************************************************
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html. Since it
//!   is executed once every logical node, it is thread-safe.
// *****************************************************************************
{
  tk::cerr_old = std::cerr.rdbuf( tk::cerr_quiet.rdbuf() );
}

QuietCerr::~QuietCerr()
// *****************************************************************************
// Destructor: restore std::cerr's stream state
// *****************************************************************************
{
  std::cerr.rdbuf( tk::cerr_old );
}

#include "NoWarning/quietcerr.def.h"
