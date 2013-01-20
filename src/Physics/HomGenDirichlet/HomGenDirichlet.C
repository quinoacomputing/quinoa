//******************************************************************************
/*!
  \file      src/Physics/HomGenDirichlet/HomGenDirichlet.C
  \author    J. Bakosi
  \date      Sat 19 Jan 2013 05:59:09 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous generalized Dirichlet model
  \details   Homogeneous generalized Dirichlet model
*/
//******************************************************************************

#include <Memory.h>
#include <MemoryException.h>
#include <MKLRandom.h>
#include <HomGenDirichlet.h>

using namespace Quinoa;

HomGenDirichlet::HomGenDirichlet(Memory* memory,
                                 Paradigm* paradigm,
                                 const int nscalar,
                                 const real time,
                                 const int echo,
                                 const int nstep) :
  Physics(memory, paradigm, "Homogeneous generalized Dirichlet", time, echo,
          nstep),
  m_nscalar(nscalar)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  nscalar  Number of mixing scalars
//! \param[in]  time     Maximum time to simulate
//! \param[in]  echo     One-line info in every few time step
//! \param[in]  nstep    Maximum number of time steps to take
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  Assert(m_random != nullptr, MemoryException,FATAL,BAD_ALLOC);
}

HomGenDirichlet::~HomGenDirichlet()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_random) { delete m_random; m_random = nullptr; }
}

void
HomGenDirichlet::echo()
//******************************************************************************
//  Echo informaion on homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
}

void
HomGenDirichlet::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
}

void
HomGenDirichlet::init()
//******************************************************************************
//  Initialize homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
}
