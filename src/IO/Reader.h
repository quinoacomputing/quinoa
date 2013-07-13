//******************************************************************************
/*!
  \file      src/IO/Reader.h
  \author    J. Bakosi
  \date      Fri 12 Jul 2013 09:33:25 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Reader base class declaration
  \details   Reader base class declaration
*/
//******************************************************************************
#ifndef Reader_h
#define Reader_h

#include <fstream>

#include <Memory.h>

namespace Quinoa {

//! Reader base
class Reader {

  protected:
    //! Constructor: Acquire file handle
    explicit Reader(const std::string filename);

    //! Destructor: Release file handle
    virtual ~Reader() noexcept;

    //! Required interface for read
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

} // namespace Quinoa

#endif // Reader_h
