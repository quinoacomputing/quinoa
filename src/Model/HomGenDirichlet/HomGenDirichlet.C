//******************************************************************************
/*!
  \file      src/Model/HomGenDirichlet/HomGenDirichlet.C
  \author    J. Bakosi
  \date      Thu Nov 15 15:05:00 2012
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
                                 const int nscalar) :
  Model(memory, paradigm, "Homogeneous generalized Dirichlet"),
  m_nscalar(nscalar)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  nscalar  Number of mixing scalars
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
  //if (m_mixModel) { delete m_mixModel; m_mixModel = nullptr; }
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
HomGenDirichlet::init()
//******************************************************************************
//  Initialize homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
  //if (m_mixModel) m_mixModel->init();
}
