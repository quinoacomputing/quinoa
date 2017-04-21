// *****************************************************************************
/*!
  \file      src/Base/ProcessControl.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     POSIX process control wrapper definitions
  \details   POSIX process control wrapper definitions.
*/
// *****************************************************************************
#ifndef ProcessControl_h
#define ProcessControl_h

#include <iosfwd>

namespace tk {

//! Remove file from file system
void rm( const std::string& file );

} // tk::

#endif // ProcessControl_h
