//******************************************************************************
/*!
  \file      src/IO/PlotWriter.C
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 08:11:50 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Plot writer base class definition
  \details   Plot writer base class definition
*/
//******************************************************************************

#include <PlotWriter.h>
#include <IOException.h>

using namespace Quinoa;

PlotWriter::PlotWriter(string filename, UnsMesh* mesh, Memory* memory) :
              m_filename(filename), m_mesh(mesh), m_memory(memory)
//******************************************************************************
//  Constructor: Acquire plot file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outPlot.open(m_filename, ofstream::out);
  if (!m_outPlot.good())
    throw IOException(FATAL, FAILED_OPEN, m_filename);
}

PlotWriter::~PlotWriter()
//******************************************************************************
//  Destructor: Release plot file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outPlot.close();
  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
  if (m_outPlot.fail())
    cerr << "WARNING: Failed to close file: " << m_filename << endl;
}
