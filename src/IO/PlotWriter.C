//******************************************************************************
/*!
  \file      src/IO/PlotWriter.C
  \author    J. Bakosi
  \date      Fri 03 May 2013 06:23:48 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Plot writer base class definition
  \details   Plot writer base class definition
*/
//******************************************************************************

#include <iostream>

#include <PlotWriter.h>
#include <Exception.h>

using namespace Quinoa;

PlotWriter::PlotWriter(const string& filename) : m_filename(filename)
//******************************************************************************
//  Constructor: Acquire plot file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outPlot.open(m_filename, ofstream::out);

  if (!m_outPlot.good())
    throw Exception(FATAL, "Failed to open file: " + m_filename);
}

PlotWriter::~PlotWriter() noexcept
//******************************************************************************
//  Destructor: Release plot file handle
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {

    m_outPlot.close();

    if (m_outPlot.fail())
      cout << "WARNING: Failed to close file: " << m_filename << endl;

  } // emit only a warning on error
    catch (exception& e) {
      cout << "WARNING: " << e.what() << endl;
    }
    catch (...) {
      cout << "UNKNOWN EXCEPTION in PlotWriter destructor" << endl
           << "Continuing anyway..." << endl;
    }
}
