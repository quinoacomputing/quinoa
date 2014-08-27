//******************************************************************************
/*!
  \file      src/Control/FileParser.C
  \author    J. Bakosi
  \date      Tue 26 Aug 2014 06:13:01 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     File parser
  \details   File parser
*/
//******************************************************************************

#include <fstream>

#include <FileParser.h>
#include <Exception.h>

using tk::FileParser;

FileParser::FileParser( const std::string& filename ) : m_filename( filename )
//******************************************************************************
//  Constructor
//! \param[in]     filename      File to parse
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream is in a valid state.
//! \author  J. Bakosi
//******************************************************************************
{
  //! Make sure there is a filename
  Assert( !filename.empty(), "No filename specified" );

  std::ifstream q;

  // Check if file exists, throw exception if it does not
  q.open( filename, std::ifstream::in );
  ErrChk( q.good(), "Failed to open file: " + filename );

  // Attempt to read a character, throw if it fails
  // It is curious that on some systems opening a directory instead of a file
  // with the above ifstream::open() call does not set the failbit. Thus we get
  // here fine, so we try to read a character from it. If it is a directory or
  // an empty file the read will fail, so we throw. Read more at: http://
  // stackoverflow.com/questions/9591036/
  // ifstream-open-doesnt-set-error-bits-when-argument-is-a-directory.
  q.get();
  ErrChk( q.good(), "Failed to read from file: " + filename );

  // Close it
  q.close();
  ErrChk( !q.fail(), "Failed to close file: " + filename );
}

void
FileParser::echoErrors( const tk::Print& print,
                        const std::vector< std::string >& errors )
//******************************************************************************
//  Echo errors accumulated during parsing
//! \param[in] errors Vector of strings of errors
//! \author  J. Bakosi
//******************************************************************************
{
//   // Underline errors
//   std::string pos( m_string.size(), ' ' );
//   for (const auto& e : errors)
//     if (e.find( "Error" ) != std::string::npos ||
//         e.find( "Warning" ) != std::string::npos) {
//       auto sloc = e.find( "at 1," );
//       if (sloc != std::string::npos) {  // if we have location info
//         // skip "at 1,"
//         sloc += 5;
//         // find a dot starting from after "at 1,"
//         const auto eloc = e.find_first_of( '.', sloc );
//         // compute location of error by extracting it from error message
//         const std::size_t errloc = std::stoi( e.substr( sloc, eloc-sloc ) ) - 1;
//         // find beginning of erroneous argument
//         sloc = m_string.rfind( ' ', errloc-1 );
//         // special-handle the first argument which has no space in front of it
//         if (sloc == std::string::npos) sloc = 0; else ++sloc;
//         // underline error
//         for (std::size_t i=sloc; i<errloc; ++i) pos[i] = '^';
//       }
//     }
// 
//   // Output errors underlined (if any) to quiet stream, list errors, and exit
//   if (!errors.empty()) {
//     print % '\n';
//     print % ">>> Parsed command line: '" % m_string % "'\n";
//     print % ">>>                       " % pos % "\n";
//     for (const auto& e : errors) print % ">>> " % e % std::endl;   // messages
//     // Exit if there were any errors
//     Throw( "Error(s) occurred, listed above, while parsing the command line\n" );
//   }
}
