//******************************************************************************
/*!
  \file      src/Model/HomDirichlet/HomDirichlet.C
  \author    J. Bakosi
  \date      Thu Nov 15 15:47:45 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Dirichlet model
  \details   Homogeneous Dirichlet model
*/
//******************************************************************************

#include <iostream>

#include <Memory.h>
#include <MemoryException.h>
#include <MKLRandom.h>
#include <HomDirichlet.h>
#include <Dirichlet.h>

using namespace Quinoa;

HomDirichlet::HomDirichlet(Memory* memory,
                           Paradigm* paradigm,
                           const int nscalar) :
  Model(memory, paradigm, "Homogeneous Dirichlet")
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

  // Instantiate Dirichlet mix model
  m_dir = new (nothrow) Dirichlet(nscalar);
  Assert(m_dir != nullptr, MemoryException,FATAL,BAD_ALLOC);
}

HomDirichlet::~HomDirichlet()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_random) { delete m_random; m_random = nullptr; }
  if (m_dir) { delete m_dir; m_dir = nullptr; }
}

void
HomDirichlet::echo()
//******************************************************************************
//  Echo informaion on homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Model: " << m_name << endl;

  // Echo information on Dirichlet mix model
  m_dir->echo();
}

void
HomDirichlet::init()
//******************************************************************************
//  Initialize homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
  //if (m_mixModel) m_mixModel->init();
}
