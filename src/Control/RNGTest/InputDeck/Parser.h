// *****************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Random number generator test suite input deck parser
  \details   This file declares the input deck, i.e., control file, parser for
    the random number generator test suite, RNGTest.
*/
// *****************************************************************************
#ifndef RNGTestInputDeckParser_h
#define RNGTestInputDeckParser_h

#include "FileParser.h"
#include "RNGTest/CmdLine/CmdLine.h"

namespace tk { class Print; }

namespace rngtest {

namespace ctr { class InputDeck; }

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
