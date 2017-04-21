// *****************************************************************************
/*!
  \file      src/Control/FileParser.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     File parser base class definition
  \details   File parser base class defintion. File parser base serves as a
    base class for various file parsers, e.g., input deck parsers. It does
    generic low-level I/O, e.g., testing whether the file to be parsed exits or
    not and associated error handling, as well as after-parser diagnostics.
*/
// *****************************************************************************

#include <map>

#include "FileParser.h"
#include "Exception.h"
#include "Reader.h"
#include "Print.h"

using tk::FileParser;

FileParser::FileParser( const std::string& filename ) : m_filename( filename )
// *****************************************************************************
//  Constructor
//! \param[in] filename File to be parsed by the parser
//! \details This constructor does basic tests in an attempt to determine if the
//!   file to be parsed exists and is in good shape and does associated error
//!   handling. This file stream is local, only used for error checking, and
//!   thus is not part of the object state here since the parser, inheriting
//!   from FileParser, e.g., walker::InputDeckParser, parses by completely
//!   outsourcing the parsing (to PEGTL), so there is no need to store the file
//!   stream handle here.
//! \author  J. Bakosi
// *****************************************************************************
{
  // Make sure there is a filename
  Assert( !filename.empty(), "No filename specified" );

  // Local file stream handle
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
FileParser::diagnostics( const tk::Print& print,
                         const std::vector< std::string >& messages )
// *****************************************************************************
//  Echo errors and warnings accumulated during parsing
//! \param[in] print    Pretty printer
//! \param[in] messages Vector of strings of errors and warnings
//! \author  J. Bakosi
// *****************************************************************************
{
  // Bundle storing multiple messages for a single errouneous line
  struct ErroneousLine {
    std::size_t dlnum;                       //!< number of digits of line num
    std::string parsed;                      //!< original line parsed
    std::string underline;                   //!< underline
    std::vector< std::string > msg;          //!< error or warning messages
    ErroneousLine() : dlnum(0), parsed(), underline(), msg() {}
  };

  Reader id( m_filename );        // file reader for extracting erroneous lines
  bool err = false;               // signaling whether there were any errors
  std::map< std::size_t, ErroneousLine > lines; // erroneous lines, key: lineno

  // Underline errors and warnings
  for (const auto& e : messages) {

    // decide if error or warning
    char underchar = ' ';
    if (e.find( "Error" ) != std::string::npos) { err = true; underchar = '^'; }
    else if (e.find( "Warning" ) != std::string::npos) underchar = '~';

    if (underchar == '^' || underchar == '~') {
      auto sloc = e.find( "at " );
      if (sloc != std::string::npos) {  // if we have location info
        // skip "at "
        sloc += 3;
        // find a comma starting from after "at "
        auto eloc = e.find_first_of( ',', sloc );
        // extract line number of error from error message
        const std::size_t lnum = std::stoul( e.substr( sloc, eloc-sloc ) );
        // store number of digits in line number
        const auto dlnum = eloc - sloc;
        // find a dot starting from after "at "
        eloc = e.find_first_of( '.', sloc );
        // skip line number
        sloc = e.find_first_of( ',', sloc ) + 1;
        // extract column number of error from error message
        const decltype(sloc) cnum = std::stoul( e.substr( sloc, eloc-sloc ) )-1;
        // store erroneous line information in map
        auto& l = lines[ lnum ];
        // store number of digits in line number
        l.dlnum = dlnum;
        // get erroneous line from file and store
        l.parsed = id.line( lnum );
        // store message
        l.msg.push_back( e );
        // start constructing underline (from scratch if first error on line)
        if (l.underline.empty()) l.underline = std::string(l.parsed.size(),' ');
        // find beginning of erroneous argument, this can be found in either e
        // (the full error message which may contain the erroneous substring
        // between single quotes and can also contain white space), or in
        // l.parsed (the erroneouss line from the file), if we find a
        // singly-quoted substring in e, we find the location of that substring
        // in l.parsed, if we don't find a singly-quoted substring in e,
        // we reverse-search l.parsed until a white space, in which case the
        // error that will be underlined will be a single word
        sloc = e.find( '\'' );
        if (sloc == std::string::npos) {
          sloc = l.parsed.rfind( ' ', cnum-1 );
        } else {
          auto sloc2 = e.find( '\'', sloc+1 );
          sloc = l.parsed.find( e.substr( sloc+1, sloc2-sloc-1 ) ) - 1;
        }
        // special-handle the beginning of the line with no space in front of it
        if (sloc == std::string::npos) sloc = 0; else ++sloc;
        // underline error and warning differently
        for (auto i=sloc; i<cnum; ++i) l.underline[i] = underchar;
      }
    }
  }

  // Output errors and warnings underlined to quiet stream and message
  for (const auto& l : lines) {
    const auto& e = l.second;
    print % '\n';
    print % ">>> Line " % l.first % ": '" % e.parsed % "'\n";
    print % ">>>" % std::string( e.dlnum+9, ' ' ) % e.underline % "\n";
    for (const auto& m : e.msg) print % ">>> " % m % '\n';
    print % '\n';
  }

  // Exit if there were any errors
  if (err) Throw( "Error(s) occurred while parsing file " + m_filename );
}
