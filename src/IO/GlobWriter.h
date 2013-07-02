//******************************************************************************
/*!
  \file      src/IO/GlobWriter.h
  \author    J. Bakosi
  \date      Tue Jul  2 15:31:35 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Glob (i.e. domain-average statistics) writer
  \details   Glob (i.e. domain-average statistics) writer
*/
//******************************************************************************
#ifndef GlobWriter_h
#define GlobWriter_h

#include <string>
#include <fstream>

#include <QuinoaTypes.h>

namespace Quinoa {

//! GlobWriter
class GlobWriter {

  public:
    //! Constructor: Acquire glob file handle
    explicit GlobWriter(std::string filename);

    //! Destructor: Release glob file handle
    ~GlobWriter() noexcept;

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

    const std::string m_filename;            //!< Glob file name
    std::ofstream m_outGlob;                 //!< Glob file output stream
};

} // namespace Quinoa

#endif // GlobWriter_h
