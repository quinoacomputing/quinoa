//******************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 07:23:19 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     UnitTest's command line parser
  \details   UnitTest's command line parser
*/
//******************************************************************************
#ifndef UnitTestCmdLineParser_h
#define UnitTestCmdLineParser_h

#include <Print.h>
#include <StringParser.h>
#include <UnitTest/CmdLine/CmdLine.h>

namespace unittest {

//! CmdLineParser : StringParser
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc, char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // unittest::

#endif // UnitTestCmdLineParser_h
