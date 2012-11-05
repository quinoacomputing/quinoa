//******************************************************************************
/*!
  \file      src/Paradigm/OpenMP.C
  \author    J. Bakosi
  \date      Sun 04 Nov 2012 09:42:07 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     OpenMP specifics
  \details   OpenMP specifics
*/
//******************************************************************************

#ifdef _OPENMP
#include <omp.h>
#endif

#include <OpenMP.h>

using namespace Quinoa;

OpenMP::OpenMP()
//******************************************************************************
//  Constructor
//! \details Query OpenMP specifics
//! \author  J. Bakosi
//******************************************************************************
{
  // Query number of OpenMP threads available
  #ifdef _OPENMP
  m_nthread = omp_get_max_threads();
  m_used = true;        // If available, OpenMP is used by default
  #else
  m_nthread = 1;
  #endif
}

bool
OpenMP::available() const
//******************************************************************************
//  Return true if compiled with OpenMP support
//! \return True if compiled with OpenMP support
//! \author  J. Bakosi
//******************************************************************************
{
  #ifdef _OPENMP
  return true;
  #else
  return false;
  #endif
}
