//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 04:00:04 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's input deck file parser
  \details   Quinoa's input deck file parser
*/
//******************************************************************************
#ifndef QuinoaInputDeckParser_h
#define QuinoaInputDeckParser_h

#include <FileParser.h>
#include <Print.h>
#include <Quinoa/CmdLine/CmdLine.h>
#include <Quinoa/InputDeck/InputDeck.h>

namespace quinoa {

//! InputDeckParser : FileParser
class InputDeckParser : public tk::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser( const tk::Print& print,
                              const ctr::CmdLine& cmdline,
                              ctr::InputDeck& inputdeck );

  private:
    //! Make requested statistics unique
    void unique( std::vector< ctr::Product >& statistics );
};

} // namespace quinoa

#endif // QuinoaInputDeckParser_h
