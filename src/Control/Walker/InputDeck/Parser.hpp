// *****************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Parser.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Walker's input deck file parser
  \details   Walker's input deck file parser
*/
// *****************************************************************************
#ifndef WalkerInputDeckParser_h
#define WalkerInputDeckParser_h

#include "FileParser.hpp"
#include "Walker/CmdLine/CmdLine.hpp"

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
