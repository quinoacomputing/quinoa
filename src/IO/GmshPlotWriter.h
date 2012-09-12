//******************************************************************************
/*!
  \file      src/IO/GmshPlotWriter.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 05:37:43 AM KST
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
    //! Don't permit copy operator
    GmshPlotWriter(const GmshPlotWriter&);
    //! Don't permit assigment operator
    GmshPlotWriter& operator=(const GmshPlotWriter&);
};

} // namespace Quinoa

#endif // GmshPlotWriter_h
