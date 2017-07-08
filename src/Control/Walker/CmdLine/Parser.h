// *****************************************************************************
/*!
  \file      src/Control/Walker/CmdLine/Parser.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Walker's command line parser
  \details   Walker's command line parser
*/
// *****************************************************************************
#ifndef WalkerCmdLineParser_h
#define WalkerCmdLineParser_h

#include "StringParser.h"
#include "Walker/CmdLine/CmdLine.h"

namespace tk { class Print; }

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
