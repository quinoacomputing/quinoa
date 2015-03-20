//******************************************************************************
/*!
  \file      src/Control/RegTest/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Fri 20 Mar 2015 12:05:49 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     RegTest's command line parser
  \details   This file declares the command-line argument parser for the
     regression test suite, RegTest.
*/
//******************************************************************************
#ifndef RegTestCmdLineParser_h
#define RegTestCmdLineParser_h

#include <Print.h>
#include <StringParser.h>
#include <RegTest/CmdLine/CmdLine.h>

namespace regtest {

//! \brief Command-line parser for RegTest.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the regression test suite, RegTest.
//! \author J. Bakosi
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc,
                            char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // regtest::

#endif // RegTestCmdLineParser_h
