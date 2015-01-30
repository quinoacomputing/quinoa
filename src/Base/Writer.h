//******************************************************************************
/*!
  \file      src/Base/Writer.h
  \author    J. Bakosi
  \date      Thu 11 Dec 2014 02:07:20 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Writer base class declaration
  \details   Writer base class declaration. Writer base serves as a base class
    for various file writers. It does generic low-level I/O, e.g., opening and
    closing a file, and associated error handling.
*/
//******************************************************************************
#ifndef Writer_h
#define Writer_h

#include <string>
#include <fstream>

#include <Exception.h>

namespace tk {

//! Writer base serves as a base class for various file writers. It does generic
//! low-level I/O, e.g., opening and closing a file, and associated error
//! handling.
class Writer {

  public:
    //! Write interface: throw if not redefined and ends up being called.
    virtual void write() { Throw( "Writer::write() is a no-op" ); }

    //! Destructor: Release file handle
    virtual ~Writer() noexcept;

  protected:
    //! Constructor: Acquire file handle. Protected: designed to be base only.
    explicit Writer( const std::string& filename,
                     std::ios_base::openmode mode = std::ios_base::out );

    const std::string m_filename;          //!< File name
    mutable std::ofstream m_outFile;       //!< File output stream
};

} // tk::

#endif // Writer_h
