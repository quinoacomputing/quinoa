//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Sat 12 Jul 2014 10:22:15 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Random number generator test suite input deck parser
  \details   Random number generator test suite input deck parser
*/
//******************************************************************************
#ifndef RNGTestInputDeckParser_h
#define RNGTestInputDeckParser_h

#include <FileParser.h>
#include <Print.h>
#include <RNGTest/CmdLine/CmdLine.h>
#include <RNGTest/InputDeck/InputDeck.h>

namespace rngtest {

//! InputDeckParser : FileParser
class InputDeckParser : public tk::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser( const tk::Print& print,
                              const ctr::CmdLine& cmdline,
                              ctr::InputDeck& inputdeck );
};

} // namespace rngtest

#endif // RNGTestInputDeckParser_h
