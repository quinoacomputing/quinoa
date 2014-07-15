//******************************************************************************
/*!
  \file      src/IO/Writer.C
  \author    J. Bakosi
  \date      Mon 14 Jul 2014 08:27:24 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Writer base class definition
  \details   Writer base class definition
*/
//******************************************************************************

#include <iostream>

#include <Writer.h>
#include <Exception.h>

using tk::Writer;

Writer::Writer(const std::string& filename) :
  m_filename(filename),
  m_outFile()
//******************************************************************************
//  Constructor: Acquire file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile.open( m_filename, std::ofstream::out );
  ErrChk( m_outFile.good(), "Failed to open file: " + m_filename );
}

Writer::~Writer() noexcept
//******************************************************************************
//  Destructor: Release file handle
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {

    m_outFile.close();

    if ( m_outFile.fail() )
      printf( ">>> WARNING: Failed to close file: %s\n", m_filename.c_str() );

  } // emit only a warning on error
    catch ( Exception& e ) {
      e.handleException();
    }
    catch ( std::exception& e ) {
      printf( ">>> WARNING: std::exception in Writer destructor: %s\n",
              e.what() );
    }
    catch (...) {
      printf( ">>> WARNING: UNKNOWN EXCEPTION in Writer destructor\n" );
    }
}
