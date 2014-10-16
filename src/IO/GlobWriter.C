//******************************************************************************
/*!
  \file      src/IO/GlobWriter.C
  \author    J. Bakosi
  \date      Wed Apr 23 11:20:25 2014
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Glob (i.e. domain-average statistics) writer
  \details   Glob (i.e. domain-average statistics) writer
*/
//******************************************************************************

#include <Macro.h>
#include <GlobWriter.h>

using quinoa::GlobWriter;

void
GlobWriter::writeGlob(const uint64_t it, const tk::real t)
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
