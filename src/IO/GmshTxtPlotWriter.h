//******************************************************************************
/*!
  \file      src/IO/GmshTxtPlotWriter.h
  \author    J. Bakosi
  \date      Tue Jul  2 15:24:20 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshTxtPlotWriter class declaration
  \details   GmshTxtPlotWriter class declaration
*/
//******************************************************************************
#ifndef GmshTxtPlotWriter_h
#define GmshTxtPlotWriter_h

#include <string>
#include <fstream>

#include <PlotWriter.h>
#include <UnsMesh.h>

namespace Quinoa {

//! GmshTxtPlotWriter : PlotWriter
class GmshTxtPlotWriter : public PlotWriter {

  public:
    //! Constructor: Acquire plot file handle
    explicit GmshTxtPlotWriter(const std::string filename,
                               UnsMesh* const mesh,
                               Memory* const memory) :
      PlotWriter(filename), m_mesh(mesh), m_memory(memory) {}

    //! Destructor: Release plot file handle
    ~GmshTxtPlotWriter() noexcept = default;

  protected:
    UnsMesh* const m_mesh;              //!< Mesh object pointer
    Memory* const m_memory;             //!< Memory object pointer

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
