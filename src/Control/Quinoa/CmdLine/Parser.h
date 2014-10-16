//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 03:52:03 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Quinoa's command line parser
  \details   Quinoa's command line parser
*/
//******************************************************************************
#ifndef QuinoaCmdLineParser_h
#define QuinoaCmdLineParser_h

#include <Print.h>
#include <StringParser.h>
#include <Quinoa/CmdLine/CmdLine.h>

namespace quinoa {

//! CmdLineParser : StringParser
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc, char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // quinoa::

#endif // QuinoaCmdLineParser_h
