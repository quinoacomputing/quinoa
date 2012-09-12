//******************************************************************************
/*!
  \file      src/IO/PlotWriter.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 04:19:20 AM KST
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

  public:
    //! Constructor: Acquire plot file handle
    PlotWriter(string filename, UnsMesh* mesh, Memory* memory);

    //! Destructor: Release plot file handle
    virtual ~PlotWriter();

  protected:
    //! Plot file name
    string m_filename;

    //! Plot file output stream
    ofstream m_outPlot;

    //! Mesh object pointer
    UnsMesh* m_mesh;

    //! Memory object pointer
    Memory* m_memory;

  private:
    //! Don't permit copy operator
    PlotWriter(const PlotWriter&);
    //! Don't permit assigment operator
    PlotWriter& operator=(const PlotWriter&);
};

} // namespace Quinoa

#endif // PlotWriter_h
