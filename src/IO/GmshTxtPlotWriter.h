//******************************************************************************
/*!
  \file      src/IO/GmshTxtPlotWriter.h
  \author    J. Bakosi
  \date      Sun 10 Mar 2013 08:16:45 PM MDT
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
class GmshTxtPlotWriter : private PlotWriter {

  public:
    //! Constructor: Acquire plot file handle
    GmshTxtPlotWriter(string filename, UnsMesh* mesh, Memory* memory) :
      PlotWriter(filename) {}

    //! Destructor: Release plot file handle
    ~GmshTxtPlotWriter() = default;

  protected:
    //! Mesh object pointer
    UnsMesh* m_mesh;

    //! Memory object pointer
    Memory* m_memory;

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
