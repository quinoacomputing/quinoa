//******************************************************************************
/*!
  \file      src/IO/GlobWriter.C
  \author    J. Bakosi
  \date      Tue May  7 12:34:06 2013
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

GlobWriter::GlobWriter(string filename) :
  m_filename(filename),
  m_outGlob()
//******************************************************************************
//  Constructor: Acquire glob file handle
//! \param[in]  filename  File name to append to
//! \author J. Bakosi
//******************************************************************************
{
  m_outGlob.open(m_filename, ofstream::out);
  ErrChk(m_outGlob.good(), FATAL, "Failed to open file: " + m_filename);
}

GlobWriter::~GlobWriter() noexcept
//******************************************************************************
//  Destructor: Release glob file handle
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
try {

  m_outGlob.close();
  ErrChk(!m_outGlob.fail(), WARNING, "Failed to close file: " + m_filename);

} // emit only a warning on error
  catch (Exception& e) {
    e.echo("WARNING");
  }
  catch (exception& e) {
    cout << ">>> std::exception in GlobWriter destructor: " << e.what()
         << endl;
  }
  catch (...) {
    cout << ">>> UNKNOWN EXCEPTION in GlobWriter destructor" << endl;
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
