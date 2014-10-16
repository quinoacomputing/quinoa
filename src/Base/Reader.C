//******************************************************************************
/*!
  \file      src/Base/Reader.C
  \author    J. Bakosi
  \date      Thu 28 Aug 2014 03:57:39 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Reader class definition
  \details   Reader class definition
*/
//******************************************************************************

#include <iostream>

#include <Reader.h>
#include <Exception.h>

using tk::Reader;

Reader::Reader( const std::string& filename ) : m_filename( filename )
//******************************************************************************
//  Constructor: Acquire file handle
//! \author J. Bakosi
//******************************************************************************
{
  //! Make sure there is a filename
  Assert( !filename.empty(), "No filename specified" );

  // Check if file exists, throw exception if it does not
  m_inFile.open( filename, std::ifstream::in );
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
//******************************************************************************
//  Destructor: Release file handle
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
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
//******************************************************************************
//  Return first line (for detection of file type based on header)
//! \author J. Bakosi
//******************************************************************************
{
  std::string s;
  std::getline( m_inFile, s );          // read the first line
  m_inFile.seekg( 0, std::ios::beg );   // seek back to the beginning of file
  return s;
}

std::vector< std::string >
Reader::lines()
//******************************************************************************
// Read file and return a string for each line
//! \author J. Bakosi
//******************************************************************************
{
  std::string s;
  std::vector< std::string > ls;
  while ( std::getline( m_inFile, s ) ) ls.emplace_back( s );
  return ls;
}

std::string
Reader::line( std::size_t lineNum )
//******************************************************************************
// Read a given line from file
//! \author J. Bakosi
//******************************************************************************
{
  std::string s;
  std::size_t num = 0;
  while ( std::getline( m_inFile, s ) && ++num < lineNum ) {}
  m_inFile.seekg( 0, std::ios::beg );   // seek back to the beginning of file
  return s;
}
