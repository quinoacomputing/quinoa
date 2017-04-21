// *****************************************************************************
/*!
  \file      src/Control/Inciter/CmdLine/Parser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Inciter's command line parser
  \details   This file declares the command-line argument parser for the
     computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************
#ifndef InciterCmdLineParser_h
#define InciterCmdLineParser_h

#include "StringParser.h"
#include "Inciter/CmdLine/CmdLine.h"

namespace tk { class Print; }

namespace inciter {

//! \brief Command-line parser for Inciter.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the computational shock hydrodynamics
//!   tool, Inciter.
//! \author J. Bakosi
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc, char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // inciter::

#endif // InciterCmdLineParser_h
