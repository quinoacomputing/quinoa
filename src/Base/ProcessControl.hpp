// *****************************************************************************
/*!
  \file      src/Base/ProcessControl.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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
