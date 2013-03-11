//******************************************************************************
/*!
  \file      src/IO/TxtPlotWriter.h
  \author    J. Bakosi
  \date      Sun 10 Mar 2013 08:50:33 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Text plot writer
  \details   Text plot writer
*/
//******************************************************************************
#ifndef TxtPlotWriter_h
#define TxtPlotWriter_h

#include <string>
#include <fstream>

#include <QuinoaTypes.h>
#include <PlotWriter.h>

using namespace std;

namespace Quinoa {

//! TxtPlotWriter : PlotWriter
class TxtPlotWriter : public PlotWriter {

  public:
    //! Constructor: Acquire plot file handle
    TxtPlotWriter(string filename);

    //! Destructor: Release plot file handle
    ~TxtPlotWriter() = default;

    //! Plot file name
    string m_filename;

    //! Plot file output stream
    ofstream m_outPlot;

    //! Write plot file
    virtual void write();

  private:
    //! Don't permit copy constructor
    TxtPlotWriter(const TxtPlotWriter&) = delete;
    //! Don't permit copy assigment
    TxtPlotWriter& operator=(const TxtPlotWriter&) = delete;
    //! Don't permit move constructor
    TxtPlotWriter(TxtPlotWriter&&) = delete;
    //! Don't permit move assigment
    TxtPlotWriter& operator=(TxtPlotWriter&&) = delete;
};

} // namespace Quinoa

#endif // TxtPlotWriter_h
