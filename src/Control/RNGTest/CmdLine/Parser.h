//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Sun 08 Jun 2014 04:05:20 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     RNGTest's command line parser
  \details   RNGTest's command line parser
*/
//******************************************************************************
#ifndef RNGTestCmdLineParser_h
#define RNGTestCmdLineParser_h

#include <Print.h>
#include <StringParser.h>
#include <RNGTest/CmdLine/CmdLine.h>

namespace rngtest {

//! CmdLineParser : StringParser
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc, char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // rngtest::

#endif // RNGTestCmdLineParser_h
