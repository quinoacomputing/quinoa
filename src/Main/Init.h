//******************************************************************************
/*!
  \file      src/Main/Init.h
  \author    J. Bakosi
  \date      Fri 04 Oct 2013 08:02:45 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Common initialization for mains
  \details   Common initialization for mains
*/
//******************************************************************************
#ifndef Init_h
#define Init_h

#include <string>

#include <Print.h>

namespace init {

//!  Wrapper for POSIX API's getcwd() from unistd.h
std::string workdir();

//!  Wrapper for the standard C library's gettimeofday() from
std::string curtime();

//!  Echo program title
void echoHeader(const quinoa::Print& print, const std::string& title);

//!  Echo build environment
void echoBuildEnv(const quinoa::Print& print, const std::string& executable);

//!  Echo runtime environment
void echoRunEnv(const quinoa::Print& print, int argc, char** argv);

} // init::

#endif // Init_h
