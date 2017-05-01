// *****************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Parser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     UnitTest's command line parser
  \details   This file declares the command-line argument parser for the unit
     test suite, UnitTest.
*/
// *****************************************************************************
#ifndef UnitTestCmdLineParser_h
#define UnitTestCmdLineParser_h

#include "StringParser.h"
#include "UnitTest/CmdLine/CmdLine.h"

namespace tk { class Print; }

namespace unittest {

//! \brief Command-line parser for UnitTest.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the unit test suite, UnitTest.
//! \author J. Bakosi
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc,
                            char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline,
                            bool& helped );
};

} // unittest::

#endif // UnitTestCmdLineParser_h
