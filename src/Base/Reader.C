// *****************************************************************************
/*!
  \file      src/Base/Reader.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Reader class definition
  \details   Reader base class declaration. Reader base servers as a base class
    for various file readers. It does generic low-level I/O, e.g., opening and
    closing a file, and associated error handling.
*/
// *****************************************************************************

#include <cstdio>
#include <exception>
#include <string>

#include "Reader.h"
#include "Exception.h"

using tk::Reader;

Reader::Reader( const std::string& filename, std::ios_base::openmode mode ) :
  m_filename( filename ), m_inFile()
// *****************************************************************************
//  Constructor: Acquire file handle
//! \param[in] filename Name of file to open for reading
//! \param[in] mode Open mode, see
//!   http://en.cppreference.com/w/cpp/io/ios_base/openmode
//! \author J. Bakosi
// *****************************************************************************
{
  // Make sure there is a filename
  Assert( !filename.empty(), "No filename specified" );

  // Check if file exists, throw exception if it does not
  m_inFile.open( filename, mode );
  ErrChk( m_inFile.good(), "Failed to open file: " + filename );

  // Attempt to read a character, throw if it fails
  // It is curious that on some systems opening a directory instead of a file
  // with the above ifstream::open() call does not set the failbit. Thus we get
  // here fine, so we try to read a character from it. If it is a directory or
  // an empty file the read will fail, so we throw. Read more at: http://
  // stackoverflow.com/questions/9591036/
  // ifstream-open-doesnt-set-error-bits-when-argument-is-a-directory.
  m_inFile.get();
  ErrChk( m_inFile.good(), "Failed to read from file: " + filename );

  // Close it
  m_inFile.close();
  ErrChk( !m_inFile.fail(), "Failed to close file: " + filename );

  // Re-open
  m_inFile.open( filename, std::ifstream::in );
  ErrChk( m_inFile.good(), "Failed to open file: " + filename );
}

Reader::~Reader() noexcept
// *****************************************************************************
//  Destructor: Release file handle
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//!   Error handling, while done by throwing and catching exceptions, results in
//!   warnings to terminal. We use C-style printf, since that will not throw
//!   exceptions.
//! \author J. Bakosi
// *****************************************************************************
{
  // Clear failbit triggered by eof, so close() won't throw a false FAILED_CLOSE
  m_inFile.clear();

  try {

    m_inFile.close();
    ErrChk( !m_inFile.fail(), "Failed to close file: " + m_filename );

  } // emit only a warning on error
    catch (Exception& e) {
      e.handleException();
    }
    catch (std::exception& e) {
      printf( ">>> WARNING: std::exception in MeshReader destructor: %s\n",
              e.what() );
    }
    catch (...) {
      printf( ">>> WARNING: UNKNOWN EXCEPTION in MeshReader destructor\n" );
    }
}

std::string
Reader::firstline()
// *****************************************************************************
//  Return first line (for detection of file type based on header)
//! \return First line read from file. This can be used for detection of file
//!   type based on header.
//! \author J. Bakosi
// *****************************************************************************
{
  std::string s;
  std::getline( m_inFile, s );          // read the first line
  m_inFile.seekg( 0, std::ios::beg );   // seek back to the beginning of file
  return s;
}

std::vector< std::string >
Reader::lines()
// *****************************************************************************
// Read file and return a string for each line
//! \return A std::vector< std::string >, a string for each line of a file.
//! \author J. Bakosi
// *****************************************************************************
{
  std::string s;
  std::vector< std::string > ls;
  while ( std::getline( m_inFile, s ) ) ls.push_back( s );
  return ls;
}

std::string
Reader::line( std::size_t lineNum )
// *****************************************************************************
// Read a given line from file
//! \param[in] lineNum Line number to read from file
//! \return Line read from file at line given
//! \author J. Bakosi
// *****************************************************************************
{
  std::string s;
  std::size_t num = 0;
  while ( std::getline( m_inFile, s ) && ++num < lineNum ) {}
  m_inFile.seekg( 0, std::ios::beg );   // seek back to the beginning of file
  return s;
}
