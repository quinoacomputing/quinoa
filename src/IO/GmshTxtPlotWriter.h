//******************************************************************************
/*!
  \file      src/IO/GmshTxtPlotWriter.h
  \author    J. Bakosi
  \date      Tue 18 Sep 2012 09:21:59 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshTxtPlotWriter class declaration
  \details   GmshTxtPlotWriter class declaration
*/
//******************************************************************************
#ifndef GmshTxtPlotWriter_h
#define GmshTxtPlotWriter_h

#include <string>
#include <fstream>

using namespace std;

#include <PlotWriter.h>
#include <UnsMesh.h>

namespace Quinoa {

//! GmshTxtPlotWriter : PlotWriter
class GmshTxtPlotWriter : PlotWriter {

  public:
    //! Constructor: Acquire plot file handle
    GmshTxtPlotWriter(string filename, UnsMesh* mesh, Memory* memory) :
      PlotWriter(filename, mesh, memory) {}

    //! Destructor: Release plot file handle
    ~GmshTxtPlotWriter() = default;

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
