//******************************************************************************
/*!
  \file      src/IO/PlotWriter.C
  \author    J. Bakosi
  \date      Fri Apr 26 17:11:49 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Plot writer base class definition
  \details   Plot writer base class definition
*/
//******************************************************************************

#include <iostream>

#include <PlotWriter.h>
#include <IOException.h>

using namespace Quinoa;

PlotWriter::PlotWriter(const string& filename) : m_filename(filename)
//******************************************************************************
//  Constructor: Acquire plot file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outPlot.open(m_filename, ofstream::out);
  Assert(m_outPlot.good(), IOException,FATAL,IO_FAILED_OPEN,m_filename);
}

PlotWriter::~PlotWriter() noexcept
//******************************************************************************
//  Destructor: Release plot file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outPlot.close();
  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
  if (m_outPlot.fail())
    cout << "WARNING: Failed to close plot file: " << m_filename << endl;
}
