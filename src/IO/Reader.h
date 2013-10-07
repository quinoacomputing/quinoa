//******************************************************************************
/*!
  \file      src/IO/Reader.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:25:16 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Reader base class declaration
  \details   Reader base class declaration
*/
//******************************************************************************
#ifndef Reader_h
#define Reader_h

#include <fstream>

namespace tk {

//! Reader base
class Reader {

  protected:
    //! Constructor: Acquire file handle
    explicit Reader(const std::string filename);

    //! Destructor: Release file handle
    virtual ~Reader() noexcept;

    //! Read interface
    virtual void read() = 0;

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
