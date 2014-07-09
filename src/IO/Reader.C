//******************************************************************************
/*!
  \file      src/IO/Reader.C
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:30:16 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Reader class definition
  \details   Reader class definition
*/
//******************************************************************************

#include <iostream>

#include <Reader.h>
#include <Exception.h>

using tk::Reader;

Reader::Reader(const std::string& filename) :
  m_filename(filename)
//******************************************************************************
//  Constructor: Acquire file handle
//! \author J. Bakosi
//******************************************************************************
{
  // Check if file exists, throw exception if it does not
  m_inFile.open(filename, std::ifstream::in);
  ErrChk (m_inFile.good(), "Failed to open file: " + filename );

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
  m_inFile.open(filename, std::ifstream::in);
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
      e.echo("WARNING");
    }
    catch (std::exception& e) {
      printf( ">>> WARNING: std::exception in MeshReader destructor: %s\n"
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
  getline( m_inFile, s );               // read the first line
  m_inFile.seekg( 0, std::ios::beg );   // seek back to the beginning of file
  return s;
}
