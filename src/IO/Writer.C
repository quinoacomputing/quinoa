//******************************************************************************
/*!
  \file      src/IO/Writer.C
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:29:29 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
      printf( ">>> WARNING: Failed to close file: %s\n", m_filename );

  } // emit only a warning on error
    catch ( Exception& e ) {
      e.echo( "WARNING" );
    }
    catch ( std::exception& e ) {
      printf( ">>> WARNING: std::exception in Writer destructor: %s\n",
              e.what() );
    }
    catch (...) {
      printf( ">>> WARNING: UNKNOWN EXCEPTION in Writer destructor\n" );
    }
}
