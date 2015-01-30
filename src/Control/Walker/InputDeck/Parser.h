//******************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 12:33:48 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Walker's input deck file parser
  \details   Walker's input deck file parser
*/
//******************************************************************************
#ifndef WalkerInputDeckParser_h
#define WalkerInputDeckParser_h

#include <FileParser.h>
#include <Print.h>
#include <Walker/CmdLine/CmdLine.h>
#include <Walker/InputDeck/InputDeck.h>

namespace walker {

//! InputDeckParser : FileParser
class InputDeckParser : public tk::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser( const tk::Print& print,
                              const ctr::CmdLine& cmdline,
                              ctr::InputDeck& inputdeck );
};

} // namespace walker

#endif // WalkerInputDeckParser_h
