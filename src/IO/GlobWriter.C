//******************************************************************************
/*!
  \file      src/IO/GlobWriter.C
  \author    J. Bakosi
  \date      Sat 04 May 2013 06:42:41 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Glob (i.e. domain-average statistics) writer
  \details   Glob (i.e. domain-average statistics) writer
*/
//******************************************************************************

#include <iostream>

#include <Macro.h>
#include <GlobWriter.h>
#include <Exception.h>

using namespace Quinoa;

GlobWriter::GlobWriter(string filename) : m_filename(filename)
//******************************************************************************
//  Constructor: Acquire glob file handle
//! \param[in]  filename  File name to append to
//! \author J. Bakosi
//******************************************************************************
{
  m_outGlob.open(m_filename, ofstream::out);

  if (!m_outGlob.good())
    throw Exception(FATAL, "Failed to open file: " + m_filename);
}

GlobWriter::~GlobWriter() noexcept
//******************************************************************************
//  Destructor: Release glob file handle
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {

    m_outGlob.close();

    if (m_outGlob.fail())
      cout << "WARNING: Failed to close file: " << m_filename << endl;

  } // emit only a warning on error
    catch (exception& e) {
      cout << "WARNING: " << e.what() << endl;
    }
    catch (...) {
      cout << "UNKNOWN EXCEPTION in GlobWriter destructor" << endl
           << "Continuing anyway..." << endl;
    }
}

void
GlobWriter::write(const int it, const real t)
//******************************************************************************
//  Write out glob file
//! \param[in]  it         Iteration counter
//! \param[in]  t          Time
//! \author J. Bakosi
//******************************************************************************
{
  IGNORE(it);
  IGNORE(t);
}
