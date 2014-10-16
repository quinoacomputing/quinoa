//******************************************************************************
/*!
  \file      src/Base/Writer.C
  \author    J. Bakosi
  \date      Thu 11 Sep 2014 11:37:50 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Writer base class definition
  \details   Writer base class definition
*/
//******************************************************************************

#include <Writer.h>
#include <Exception.h>

using tk::Writer;

Writer::Writer( const std::string& filename, std::ios_base::openmode mode ) :
  m_filename( filename )
//******************************************************************************
//  Constructor: Acquire file handle
//! \author J. Bakosi
//******************************************************************************
{
  // Doing this if-test gives the derived class an option to pass an empty
  // string in case the file does not need to be opened, because e.g., there is
  // no data that needs to be written, without contaminating client-code with
  // this if-test.
  if (!m_filename.empty()) {
    m_outFile.open( m_filename, mode );
    ErrChk( m_outFile.good(), "Failed to open file: " + m_filename );
  }
}

Writer::~Writer() noexcept
//******************************************************************************
//  Destructor: Release file handle
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  if (!m_filename.empty()) {
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
}
