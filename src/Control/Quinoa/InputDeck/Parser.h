//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 12:33:51 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Quinoa's input deck file parser
  \details   This file declares the input deck, i.e., control file, parser for
    the computational fluid dynamics tool, Quinoa.
*/
//******************************************************************************
#ifndef QuinoaInputDeckParser_h
#define QuinoaInputDeckParser_h

#include <FileParser.h>
#include <Print.h>
#include <Quinoa/CmdLine/CmdLine.h>
#include <Quinoa/InputDeck/InputDeck.h>

namespace quinoa {

//! \brief Control file parser for Quinoa.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing the control file for the computational fluid dynamics tool,
//!   Quinoa.
class InputDeckParser : public tk::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser( const tk::Print& print,
                              const ctr::CmdLine& cmdline,
                              ctr::InputDeck& inputdeck );
};

} // namespace quinoa

#endif // QuinoaInputDeckParser_h
