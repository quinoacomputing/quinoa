//******************************************************************************
/*!
  \file      HomogeneousGeneralizedDirichlet.C
  \author    J. Bakosi
  \date      Sun 20 Jan 2013 01:54:26 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous generalized Dirichlet model
  \details   Homogeneous generalized Dirichlet model
*/
//******************************************************************************

#include <Memory.h>
#include <MemoryException.h>
#include <MKLRandom.h>
#include <HomogeneousGeneralizedDirichlet.h>

using namespace Quinoa;

HomogeneousGeneralizedDirichlet::HomogeneousGeneralizedDirichlet(
                                   Memory* memory,
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

HomogeneousGeneralizedDirichlet::~HomogeneousGeneralizedDirichlet()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_random) { delete m_random; m_random = nullptr; }
}

void
HomogeneousGeneralizedDirichlet::echo()
//******************************************************************************
//  Echo informaion on homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
}

void
HomogeneousGeneralizedDirichlet::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
}

void
HomogeneousGeneralizedDirichlet::init()
//******************************************************************************
//  Initialize homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
}
