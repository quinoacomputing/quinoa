//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Keywords.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 09:35:20 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's command line keywords
  \details   All keywords recognized by Quinoa's command line parser. The
  keywords are defined by specializing struct 'keyword', defined in
  Control/Keyword.h. Introducing a new keyword requires a more human readable
  (but still short) name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef QuinoaCmdLineKeywords_h
#define QuinoaCmdLineKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes below to make sure they get included in the correct
//! namespace and not polluting the global one.
#define Keywords



#undef Keywords

#endif // QuinoaCmdLineKeywords_h
