//******************************************************************************
/*!
  \file      src/Base/Reader.h
  \author    J. Bakosi
  \date      Thu 28 Aug 2014 12:20:07 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Reader base class declaration
  \details   Reader base class declaration
*/
//******************************************************************************
#ifndef Reader_h
#define Reader_h

#include <fstream>
#include <vector>

#include <Exception.h>

namespace tk {

//! Reader base
class Reader {

  public:
    //! Constructor: Acquire file handle
    explicit Reader( const std::string& filename );

    //! Destructor: Release file handle
    virtual ~Reader() noexcept;

    //! Read interface
    virtual void read() { Throw( "Reader::read() is a no-op" ); }

    //! Return first line (for detection of file type based on header)
    std::string firstline();

    //! Read a given line from file
    std::string line( std::size_t lineNum );

    //! Read file and return a string for each line
    std::vector< std::string > lines();

  protected:
    const std::string m_filename;            //!< File name
    std::ifstream m_inFile;                  //!< File input stream
};

} // tk::

#endif // Reader_h
