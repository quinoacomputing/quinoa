//******************************************************************************
/*!
  \file      src/Base/Reader.h
  \author    J. Bakosi
  \date      Wed 18 Feb 2015 09:54:02 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Reader base class declaration
  \details   Reader base class declaration. Reader base serves as a base class
    for various file readers. It does generic low-level I/O, e.g., opening and
    closing a file, and associated error handling.
*/
//******************************************************************************
#ifndef Reader_h
#define Reader_h

#include <fstream>
#include <vector>

#include <Exception.h>

namespace tk {

//! Reader base serves as a base class for various file readers. It does generic
//! low-level I/O, e.g., opening and closing a file, and associated error
//! handling.
class Reader {

  public:
    //! Constructor: Acquire file handle
    explicit Reader( const std::string& filename );

    //! Destructor: Release file handle
    virtual ~Reader() noexcept;

    //! Read interface: throw if not redefined and ends up being called called.
    virtual void read() { Throw( "Reader::read() is a no-op" ); }

    //! Return first line (for detection of file type based on header)
    std::string firstline();

    //! Read file and return a string for each line
    std::vector< std::string > lines();

    //! Read a given line from file
    std::string line( std::size_t lineNum );

  protected:
    const std::string m_filename;            //!< File name
    std::ifstream m_inFile;                  //!< File input stream
};

} // tk::

#endif // Reader_h
