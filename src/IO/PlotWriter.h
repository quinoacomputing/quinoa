//******************************************************************************
/*!
  \file      src/IO/PlotWriter.h
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 05:59:03 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     PlotWriter base class declaration
  \details   PlotWriter base class declaration
*/
//******************************************************************************
#ifndef PlotWriter_h
#define PlotWriter_h

#include <string>
#include <fstream>

using namespace std;

#include <UnsMesh.h>

namespace Quinoa {

//! PlotWriter base class
class PlotWriter {

  protected:
    //! Constructor: Acquire plot file handle
    PlotWriter(string filename, UnsMesh* mesh, Memory* memory);

    //! Destructor: Release plot file handle
    virtual ~PlotWriter();

    //! Plot file name
    string m_filename;

    //! Plot file output stream
    ofstream m_outPlot;

    //! Mesh object pointer
    UnsMesh* m_mesh;

    //! Memory object pointer
    Memory* m_memory;

  private:
    //! Don't permit copy constructor
    PlotWriter(const PlotWriter&) = delete;
    //! Don't permit copy assigment
    PlotWriter& operator=(const PlotWriter&) = delete;
    //! Don't permit move constructor
    PlotWriter(PlotWriter&&) = delete;
    //! Don't permit move assigment
    PlotWriter& operator=(PlotWriter&&) = delete;
};

} // namespace Quinoa

#endif // PlotWriter_h
