//******************************************************************************
/*!
  \file      src/IO/GlobWriter.h
  \author    J. Bakosi
  \date      Sun 10 Mar 2013 08:44:37 PM MDT
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

using namespace std;

namespace Quinoa {

//! GlobWriter
class GlobWriter {

  public:
    //! Constructor: Acquire glob file handle
    GlobWriter(string filename);

    //! Destructor: Release glob file handle
    ~GlobWriter();

    //! Write glob file
    void write();

  private:
    //! Don't permit copy constructor
    GlobWriter(const GlobWriter&) = delete;
    //! Don't permit copy assigment
    GlobWriter& operator=(const GlobWriter&) = delete;
    //! Don't permit move constructor
    GlobWriter(GlobWriter&&) = delete;
    //! Don't permit move assigment
    GlobWriter& operator=(GlobWriter&&) = delete;

    const string m_filename;            //!< Glob file name
    ofstream m_outGlob;                 //!< Glob file output stream
};

} // namespace Quinoa

#endif // GlobWriter_h
