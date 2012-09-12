//******************************************************************************
/*!
  \file      src/IO/PlotWriter.C
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 04:42:14 AM KST
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
  if (!m_outPlot.good()) throw IOException(FATAL, FAILED_OPEN, m_filename);
}

PlotWriter::~PlotWriter()
//******************************************************************************
//  Destructor: Release plot file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outPlot.close();
  if (m_outPlot.fail()) throw IOException(WARNING, FAILED_CLOSE, m_filename);
}
