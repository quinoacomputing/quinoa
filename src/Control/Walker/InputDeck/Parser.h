// *****************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Parser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Walker's input deck file parser
  \details   Walker's input deck file parser
*/
// *****************************************************************************
#ifndef WalkerInputDeckParser_h
#define WalkerInputDeckParser_h

#include "FileParser.h"
#include "Walker/CmdLine/CmdLine.h"

namespace tk { class Print; }

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
