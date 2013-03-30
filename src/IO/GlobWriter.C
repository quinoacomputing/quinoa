//******************************************************************************
/*!
  \file      src/IO/GlobWriter.C
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 11:16:22 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Glob (i.e. domain-average statistics) writer
  \details   Glob (i.e. domain-average statistics) writer
*/
//******************************************************************************

#include <iostream>

#include <Macro.h>
#include <GlobWriter.h>
#include <IOException.h>

using namespace Quinoa;

GlobWriter::GlobWriter(string filename) : m_filename(filename)
//******************************************************************************
//  Constructor: Acquire glob file handle
//! \param[in]  filename  File name to append to
//! \author J. Bakosi
//******************************************************************************
{
  m_outGlob.open(m_filename, ofstream::out);
  Assert(m_outGlob.good(), IOException,FATAL,IO_FAILED_OPEN,m_filename);
}

GlobWriter::~GlobWriter()
//******************************************************************************
//  Destructor: Release glob file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outGlob.close();
  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
  if (m_outGlob.fail())
    cout << "WARNING: Failed to close glob file: " << m_filename << endl;
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
