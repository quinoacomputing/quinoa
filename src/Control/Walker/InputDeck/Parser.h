//******************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 06:36:52 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
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

  private:
    //! Make requested statistics unique
    void unique( std::vector< tk::ctr::Product >& statistics );
};

} // namespace walker

#endif // WalkerInputDeckParser_h
