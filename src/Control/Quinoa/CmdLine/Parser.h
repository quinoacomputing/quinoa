//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 06:17:13 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Quinoa's command line parser
  \details   This file declares the command-line argument parser for the
     computational fluid dynamics tool, Quinoa.
*/
//******************************************************************************
#ifndef QuinoaCmdLineParser_h
#define QuinoaCmdLineParser_h

#include <Print.h>
#include <StringParser.h>
#include <Quinoa/CmdLine/CmdLine.h>

namespace quinoa {

//! \brief Command-line parser for Quinoa.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the computational fluid dynamics tool,
//!   Quinoa.
//! \author J. Bakosi
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc,
                            char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // quinoa::

#endif // QuinoaCmdLineParser_h
