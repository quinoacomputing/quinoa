//******************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Sun 12 Apr 2015 07:22:44 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter's input deck file parser
  \details   This file declares the input deck, i.e., control file, parser for
    the computational shock hydrodynamics tool, Inciter.
*/
//******************************************************************************
#ifndef InciterInputDeckParser_h
#define InciterInputDeckParser_h

#include <FileParser.h>
#include <Print.h>
#include <Inciter/CmdLine/CmdLine.h>
#include <Inciter/InputDeck/InputDeck.h>

namespace inciter {

//! \brief Control file parser for Inciter.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing the control file for the computational shock hydrodynamics tool,
//!   Inciter.
class InputDeckParser : public tk::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser( const tk::Print& print,
                              const ctr::CmdLine& cmdline,
                              ctr::InputDeck& inputdeck,
                              int peid );
};

} // namespace inciter

#endif // InciterInputDeckParser_h
