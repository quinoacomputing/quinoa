//******************************************************************************
/*!
  \file      src/IO/GlobWriter.C
  \author    J. Bakosi
  \date      Thu Aug 29 15:32:54 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Glob (i.e. domain-average statistics) writer
  \details   Glob (i.e. domain-average statistics) writer
*/
//******************************************************************************

#include <Macro.h>
#include <GlobWriter.h>

using namespace quinoa;

void
GlobWriter::write(const uint64_t it, const real t)
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
