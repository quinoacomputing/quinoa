//******************************************************************************
/*!
  \file      src/Main/Init.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:30:46 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Common initialization for mains
  \details   Common initialization for mains
*/
//******************************************************************************
#ifndef Init_h
#define Init_h

#include <string>

#include <Print.h>

namespace tk {

//!  Wrapper for POSIX API's getcwd() from unistd.h
std::string workdir();

//!  Wrapper for the standard C library's gettimeofday() from
std::string curtime();

//!  Echo program title
void echoHeader(const tk::Print& print, const std::string& title);

//!  Echo build environment
void echoBuildEnv(const tk::Print& print, const std::string& executable);

//!  Echo runtime environment
void echoRunEnv(const tk::Print& print, int argc, char** argv);

} // tk::

#endif // Init_h
