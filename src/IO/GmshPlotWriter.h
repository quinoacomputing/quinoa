//******************************************************************************
/*!
  \file      src/IO/GmshPlotWriter.h
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 05:22:29 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshPlotWriter class declaration
  \details   GmshPlotWriter class declaration
*/
//******************************************************************************
#ifndef GmshPlotWriter_h
#define GmshPlotWriter_h

#include <string>
#include <fstream>

using namespace std;

#include <PlotWriter.h>
#include <UnsMesh.h>

namespace Quinoa {

//! GmshPlotWriter base class
class GmshPlotWriter : PlotWriter {

  public:
    //! Constructor: Acquire plot file handle
    GmshPlotWriter(string filename, UnsMesh* mesh, Memory* memory) :
      PlotWriter(filename, mesh, memory) {}

    //! Destructor: Release plot file handle
    virtual ~GmshPlotWriter();

  private:
    //! Don't permit copy constructor
    GmshPlotWriter(const GmshPlotWriter&) = delete;
    //! Don't permit copy assigment
    GmshPlotWriter& operator=(const GmshPlotWriter&) = delete;
    //! Don't permit move constructor
    GmshPlotWriter(GmshPlotWriter&&) = delete;
    //! Don't permit move assigment
    GmshPlotWriter& operator=(GmshPlotWriter&&) = delete;
};

} // namespace Quinoa

#endif // GmshPlotWriter_h
