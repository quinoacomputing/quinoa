//******************************************************************************
/*!
  \file      src/Control/RegTest/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Fri 29 May 2015 11:41:35 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     RegTest's command line parser
  \details   This file declares the command-line argument parser for the
     regression test suite, RegTest.
*/
//******************************************************************************
#ifndef RegTestCmdLineParser_h
#define RegTestCmdLineParser_h

#include "StringParser.h"
#include "RegTest/CmdLine/CmdLine.h"

namespace tk { class Print; }

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
