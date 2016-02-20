//******************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Sat 30 May 2015 12:03:16 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Walker's input deck file parser
  \details   Walker's input deck file parser
*/
//******************************************************************************
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
