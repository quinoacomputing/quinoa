//******************************************************************************
/*!
  \file      src/Control/StringParser.C
  \author    J. Bakosi
  \date      Sun 10 Nov 2013 06:42:01 AM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     String parser
  \details   String parser
*/
//******************************************************************************

#include <StringParser.h>

using tk::StringParser;

StringParser::StringParser(int argc, char** argv)
//******************************************************************************
//  Constructor
//! \param[in]     argc          Number of C-style character arrays in argv
//! \param[in]     argv          C-style character array of character arrays
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream is in a valid state.
//! \author  J. Bakosi
//******************************************************************************
{
  // Convert C-style character array of character arrays to a single std::string
  // substrings separated by spaces
  for (int i=1; i<argc; ++i) {
    m_string += std::string(argv[i]) + ' ';
  }
}
