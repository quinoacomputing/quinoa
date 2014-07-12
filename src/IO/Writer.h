//******************************************************************************
/*!
  \file      src/IO/Writer.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 08:51:23 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Writer base class declaration
  \details   Writer base class declaration
*/
//******************************************************************************
#ifndef Writer_h
#define Writer_h

#include <string>
#include <fstream>

#include <Exception.h>

namespace tk {

//! Writer base class
class Writer {

  public:
    //! Constructor: Acquire file handle
    explicit Writer(const std::string& filename);

    //! Destructor: Release file handle
    virtual ~Writer() noexcept;

    //! Write interface
    virtual void write() {
      Throw( "Writer::write() is a no-op and should be defined by a child" );
    }

  protected:
    const std::string m_filename;          //!< File name

    mutable std::ofstream m_outFile;       //!< File output stream

  private:
    //! Don't permit copy constructor
    Writer(const Writer&) = delete;
    //! Don't permit copy assigment
    Writer& operator=(const Writer&) = delete;
    //! Don't permit move constructor
    Writer(Writer&&) = delete;
    //! Don't permit move assigment
    Writer& operator=(Writer&&) = delete;
};

} // tk::

#endif // Writer_h
