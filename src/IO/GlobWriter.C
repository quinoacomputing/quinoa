//******************************************************************************
/*!
  \file      src/IO/GlobWriter.C
  \author    J. Bakosi
  \date      Sun 10 Mar 2013 04:31:09 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Glob (i.e. domain-average statistics) writer
  \details   Glob (i.e. domain-average statistics) writer
*/
//******************************************************************************

#include <iostream>
#include <iomanip>

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

void
GlobWriter::write(const int it,
                  const real t,
                  const int nord,
                  const real* const ordinary)
//******************************************************************************
//  Write out domain-averaged statistics to file
//! \param[in]  it         Iteration counter
//! \param[in]  t          Time
//! \param[in]  time       Time stamp
//! \param[in]  nord       Number of ordinary moments
//! \param[in]  ordinary   Ordinary moments
//! \author  J. Bakosi
//******************************************************************************
{
  m_outGlob << setfill(' ') << setw(8) << it << "  " << scientific
            << setprecision(6) << setw(12) << t << "  ";

  for (int i=0; i<nord; ++i) {
    m_outGlob << ordinary[i] << "  ";
  }

  m_outGlob << endl;
}
