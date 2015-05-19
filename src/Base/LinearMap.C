//******************************************************************************
/*!
  \file      src/Base/LinearMap.C
  \author    J. Bakosi
  \date      Sat 04 Apr 2015 08:06:28 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Advanced Charm++ array creation with a map in a linear fashion
  \details   Advanced Charm++ array creation refers to various ways arrays can
     be created with the Charm++ runtime system. See
     http://charm.cs.illinois.edu/manuals/html/charm++/manual.html, Sec.
     Advanced Array Creation. This class does a simple linear distribution.
*/
//******************************************************************************

#include <LinearMap.h>
#include <Exception.h>

using tk::LinearMap;

int
LinearMap::procNum( int, const CkArrayIndex& idx )
//******************************************************************************
//  Return the home processor number for the array element for linear
//  distribution
//! \param[in] unused
//! \param[in] idx Charm++ array index object containing the array element index
//!   to assign a PE to
//! \return PE assigned
//! \author J. Bakosi
//******************************************************************************
{
  int elem = *idx.data();
  auto pe = elem / m_chunksize;
  if (pe >= CkNumPes()) pe = CkNumPes()-1;

  Assert( pe < CkNumPes(), "Assigned PE (" + std::to_string(pe) +
          ") larger than NumPEs (" + std::to_string(CkNumPes()) + ")" );

  return pe;
}

void
LinearMap::populateInitial( int, CkArrayIndex& idx, void *msg, CkArrMgr *mgr )
//******************************************************************************
// Create initial set of array elements based on linear distribution
//! \param[in] unused
//! \param[in] idx Charm++ array index object containing the number of initial
//!   array elements to be created
//! \author J. Bakosi
//******************************************************************************
{
  int nelem = *idx.data();
  if (nelem == 0) return;     // no initial elements requested

  auto lower = CkMyPe() * m_chunksize;
  auto upper = lower + m_chunksize;
  auto remainder = nelem % CkNumPes();
  if (remainder && CkMyPe() == CkNumPes()-1) upper += remainder;

  for (int e=0; e<nelem; ++e)
    if (e >= lower && e < upper)
      mgr->insertInitial( e, CkCopyMsg(&msg) );

  mgr->doneInserting();
  CkFreeMsg( msg );
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <linearmap.def.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
