//******************************************************************************
/*!
  \file      src/Paradigm/OpenMP.C
  \author    J. Bakosi
  \date      Wed 01 May 2013 09:04:13 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     OpenMP specifics
  \details   OpenMP specifics
*/
//******************************************************************************

#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <OpenMP.h>

using namespace Quinoa;

OpenMP::OpenMP() :
#ifdef _OPENMP
  m_available(true),
  m_used(true),
  m_nthread(omp_get_max_threads())
#else  // _OPENMP
  m_available(false),
  m_used(false),
  m_nthread(1)
#endif // _OPENMP
//******************************************************************************
//  Constructor
//! \details Query OpenMP specifics
//! \author  J. Bakosi
//******************************************************************************
{
  assert(m_nthread > 0);
}
