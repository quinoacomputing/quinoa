//******************************************************************************
/*!
  \file      src/IO/Reader.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:02:06 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Reader base class declaration
  \details   Reader base class declaration
*/
//******************************************************************************
#ifndef Reader_h
#define Reader_h

#include <fstream>

#include <Exception.h>

namespace tk {

//! Reader base
class Reader {

  public:
    //! Return first line (for detection of file type based on header)
    std::string firstline();

    //! Constructor: Acquire file handle
    explicit Reader(const std::string& filename);

    //! Destructor: Release file handle
    virtual ~Reader() noexcept;

    //! Read interface
    virtual void read() { Throw( "Reader::read() is a no-op"); }

  protected:
    const std::string m_filename;            //!< File name
    std::ifstream m_inFile;                  //!< File input stream

  private:
    //! Don't permit copy constructor
    Reader(const Reader&) = delete;
    //! Don't permit copy assigment
    Reader& operator=(const Reader&) = delete;
    //! Don't permit move constructor
    Reader(Reader&&) = delete;
    //! Don't permit move assigment
    Reader& operator=(Reader&&) = delete;
};

} // tk::

#endif // Reader_h
