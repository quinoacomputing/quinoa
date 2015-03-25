//******************************************************************************
/*!
  \file      src/Base/ProcessException.h
  \author    J. Bakosi
  \date      Tue 24 Mar 2015 12:23:53 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Process an exception definition
  \details   This file contains the definition of processing an exception.
    Logically, it would make sense to put this into Exception.C, however,
    Exception.h is included by all who want to be able throw an exception (a
    potentially large number of files) and that would pull in the charm++.h
    header, which triggers a slew of compiler warnings. On the other hand,
    processing an exception is only done by executable's main chares objects (a
    lot less than those which throw), so processing an exception is separated
    here.
*/
//******************************************************************************
#ifndef ProcessException_h
#define ProcessException_h

namespace tk {

//! Process an exception
void processException();

} // tk::

#endif // ProcessException_h
