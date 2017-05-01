// *****************************************************************************
/*!
  \file      src/Base/Writer.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Writer base class declaration
  \details   Writer base class declaration. Writer base serves as a base class
    for various file writers. It does generic low-level I/O, e.g., opening and
    closing a file, and associated error handling.
*/
// *****************************************************************************
#ifndef Writer_h
#define Writer_h

#include <fstream>

#include "Exception.h"

namespace tk {

//! Writer base serves as a base class for various file writers. It does generic
//! low-level I/O, e.g., opening and closing a file, and associated error
//! handling.
class Writer {

  public:
    //! Constructor: Acquire file handle. Protected: designed to be base only.
    explicit Writer( const std::string& filename,
                     std::ios_base::openmode mode = std::ios_base::out );

    //! Destructor: Release file handle
    virtual ~Writer() noexcept;

    //! Unformatted write
    //! \param[in] data Buffer to write
    //! \param[in] count Number of characters to write
    void write( const char* data, std::streamsize count )
    { m_outFile.write( data, count ); }

  protected:
    const std::string m_filename;          //!< File name
    mutable std::ofstream m_outFile;       //!< File output stream
};

} // tk::

#endif // Writer_h
