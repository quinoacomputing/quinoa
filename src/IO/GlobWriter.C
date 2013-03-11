//******************************************************************************
/*!
  \file      src/IO/GlobWriter.C
  \author    J. Bakosi
  \date      Sun 10 Mar 2013 08:47:08 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Glob (i.e. domain-average statistics) writer
  \details   Glob (i.e. domain-average statistics) writer
*/
//******************************************************************************

#include <iostream>

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
    cout << "WARNING: Failed to close file: " << m_filename << endl;
}
