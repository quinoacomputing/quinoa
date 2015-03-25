//******************************************************************************
/*!
  \file      src/Base/ProcessControl.h
  \author    J. Bakosi
  \date      Tue 24 Mar 2015 12:32:57 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     POSIX process control wrapper definitions
  \details   POSIX process control wrapper definitions.
*/
//******************************************************************************
#ifndef ProcessControl_h
#define ProcessControl_h

namespace tk {

//! Remove file from file system
void rm( const std::string& file );

} // tk::

#endif // ProcessControl_h
