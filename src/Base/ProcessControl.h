// *****************************************************************************
/*!
  \file      src/Base/ProcessControl.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
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
