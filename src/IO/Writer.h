//******************************************************************************
/*!
  \file      src/IO/Writer.h
  \author    J. Bakosi
  \date      Tue 31 Dec 2013 01:01:09 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Writer base class declaration
  \details   Writer base class declaration
*/
//******************************************************************************
#ifndef Writer_h
#define Writer_h

#include <string>
#include <fstream>

namespace tk {

//! Writer base class
class Writer {

  protected:
    //! Constructor: Acquire file handle
    explicit Writer(const std::string& filename);

    //! Destructor: Release file handle
    virtual ~Writer() noexcept;

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
