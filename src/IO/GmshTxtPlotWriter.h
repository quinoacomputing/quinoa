//******************************************************************************
/*!
  \file      src/IO/GmshTxtPlotWriter.h
  \author    J. Bakosi
  \date      Fri Jul 19 15:57:13 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshTxtPlotWriter class declaration
  \details   GmshTxtPlotWriter class declaration
*/
//******************************************************************************
#ifndef GmshTxtPlotWriter_h
#define GmshTxtPlotWriter_h

#include <string>

#include <Writer.h>

namespace Quinoa {

//! GmshTxtPlotWriter : Writer
class GmshTxtPlotWriter : public Writer {

  public:
    //! Constructor: Acquire plot file handle
    explicit GmshTxtPlotWriter(const std::string& filename) :
      Writer(filename) {}

    //! Destructor: Release plot file handle
    ~GmshTxtPlotWriter() noexcept = default;

  private:
    //! Don't permit copy constructor
    GmshTxtPlotWriter(const GmshTxtPlotWriter&) = delete;
    //! Don't permit copy assigment
    GmshTxtPlotWriter& operator=(const GmshTxtPlotWriter&) = delete;
    //! Don't permit move constructor
    GmshTxtPlotWriter(GmshTxtPlotWriter&&) = delete;
    //! Don't permit move assigment
    GmshTxtPlotWriter& operator=(GmshTxtPlotWriter&&) = delete;
};

} // namespace Quinoa

#endif // GmshTxtPlotWriter_h
