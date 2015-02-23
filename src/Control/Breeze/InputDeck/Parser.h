//******************************************************************************
/*!
  \file      src/Control/Breeze/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 08:48:13 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Breeze's input deck file parser
  \details   This file declares the input deck, i.e., control file, parser for
    the computational fluid dynamics tool, Breeze.
*/
//******************************************************************************
#ifndef BreezeInputDeckParser_h
#define BreezeInputDeckParser_h

#include <FileParser.h>
#include <Print.h>
#include <Breeze/CmdLine/CmdLine.h>
#include <Breeze/InputDeck/InputDeck.h>

namespace inciter {

//! \brief Control file parser for Breeze.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing the control file for the computational fluid dynamics tool,
//!   Breeze.
class InputDeckParser : public tk::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser( const tk::Print& print,
                              const ctr::CmdLine& cmdline,
                              ctr::InputDeck& inputdeck );
};

} // namespace inciter

#endif // BreezeInputDeckParser_h
