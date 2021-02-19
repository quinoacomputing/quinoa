// *****************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Random number generator test suite input deck parser
  \details   This file declares the input deck, i.e., control file, parser for
    the random number generator test suite, RNGTest.
*/
// *****************************************************************************
#ifndef RNGTestInputDeckParser_h
#define RNGTestInputDeckParser_h

#include "FileParser.hpp"
#include "RNGTest/CmdLine/CmdLine.hpp"

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
