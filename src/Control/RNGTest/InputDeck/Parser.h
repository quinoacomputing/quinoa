//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 06:32:26 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Random number generator test suite input deck parser
  \details   This file declares the input deck, i.e., control file, parser for
    the random number generator test suite, RNGTest.
*/
//******************************************************************************
#ifndef RNGTestInputDeckParser_h
#define RNGTestInputDeckParser_h

#include <FileParser.h>
#include <Print.h>
#include <RNGTest/CmdLine/CmdLine.h>
#include <RNGTest/InputDeck/InputDeck.h>

namespace rngtest {

//! \brief Control file parser for RNGTest.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing the control file for the random number generator test suite,
//!   RNGTest.
class InputDeckParser : public tk::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser( const tk::Print& print,
                              const ctr::CmdLine& cmdline,
                              ctr::InputDeck& inputdeck );
};

} // namespace rngtest

#endif // RNGTestInputDeckParser_h
