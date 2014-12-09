//******************************************************************************
/*!
  \file      src/Control/Walker/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 11:31:10 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Walker's command line parser
  \details   Walker's command line parser
*/
//******************************************************************************
#ifndef WalkerCmdLineParser_h
#define WalkerCmdLineParser_h

#include <Print.h>
#include <StringParser.h>
#include <Walker/CmdLine/CmdLine.h>

namespace walker {

//! CmdLineParser : StringParser
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc, char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // walker::

#endif // WalkerCmdLineParser_h
