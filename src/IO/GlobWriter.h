//******************************************************************************
/*!
  \file      src/IO/GlobWriter.h
  \author    J. Bakosi
  \date      Fri 12 Jul 2013 10:01:22 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Glob (i.e. domain-average statistics) writer
  \details   Glob (i.e. domain-average statistics) writer
*/
//******************************************************************************
#ifndef GlobWriter_h
#define GlobWriter_h

#include <string>

#include <QuinoaTypes.h>
#include <Writer.h>

namespace Quinoa {

//! GlobWriter : Writer
class GlobWriter : public Writer {

  public:
    //! Constructor: Acquire glob file handle
    explicit GlobWriter(std::string filename) :
      Writer(filename) {}

    //! Destructor: Release glob file handle
    virtual ~GlobWriter() noexcept = default;

    //! Write glob file
    void write(const uint64_t it, const real t);

  private:
    //! Don't permit copy constructor
    GlobWriter(const GlobWriter&) = delete;
    //! Don't permit copy assigment
    GlobWriter& operator=(const GlobWriter&) = delete;
    //! Don't permit move constructor
    GlobWriter(GlobWriter&&) = delete;
    //! Don't permit move assigment
    GlobWriter& operator=(GlobWriter&&) = delete;
};

} // namespace Quinoa

#endif // GlobWriter_h
