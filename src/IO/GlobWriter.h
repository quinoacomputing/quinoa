//******************************************************************************
/*!
  \file      src/IO/GlobWriter.h
  \author    J. Bakosi
  \date      Wed Apr 23 11:20:15 2014
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Glob (i.e. domain-average statistics) writer
  \details   Glob (i.e. domain-average statistics) writer
*/
//******************************************************************************
#ifndef GlobWriter_h
#define GlobWriter_h

#include <string>

#include <Types.h>
#include <Writer.h>

namespace quinoa {

//! GlobWriter : Writer
class GlobWriter : public tk::Writer {

  public:
    //! Constructor: Acquire glob file handle
    explicit GlobWriter(const std::string& filename) :
      Writer(filename) {}

    //! Destructor: Release glob file handle
    ~GlobWriter() noexcept override = default;

    //! Write glob file
    void writeGlob(const uint64_t it, const tk::real t);

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

} // quinoa::

#endif // GlobWriter_h
