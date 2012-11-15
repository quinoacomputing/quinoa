//******************************************************************************
/*!
  \file      src/Model/Model.C
  \author    J. Bakosi
  \date      Thu Nov 15 13:38:16 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************

#include <Model.h>
#include <ModelException.h>
#include <MemoryException.h>
#include <Dirichlet.h>
#include <GeneralizedDirichlet.h>
#include <Setup.h>
#include <Memory.h>
#include <MKLRandom.h>

using namespace Quinoa;

Model::Model(const ModelType model,
             const int npel,
             Memory* memory,
             Paradigm* paradigm) :
  m_model(model), m_npel(npel), m_memory(memory), m_paradigm(paradigm)
//******************************************************************************
//  Constructor
//! \param[in]  model    Model type (see Control/Control.h)
//! \param[in]  npel     Number of particles/element
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate random number generator (the same object can then be used by
  // all models)
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  Assert(m_random != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Instantiate model
  // ICC: this could be a switch
  if (model == ModelType::HOMOGENEOUS_DIRICHLET) {
    m_name = "Homogeneous Dirichlet";
    m_nel = 1;
    m_mixModel = new (nothrow) Dirichlet(this, m_random, NSCALAR);
    Assert(m_mixModel != nullptr, MemoryException,FATAL,BAD_ALLOC);
  }
  else if (model == ModelType::HOMOGENEOUS_GENDIRICHLET) {
    m_name = "Homogeneous generalized Dirichlet";
    m_nel = 1;
    m_mixModel = new (nothrow) GeneralizedDirichlet(this, m_random, NSCALAR);
    Assert(m_mixModel != nullptr, MemoryException,FATAL,BAD_ALLOC);
  } else {
    Throw(ModelException,FATAL,NO_SUCH_MODEL);
  }

  // Zero out MemoryEntry pointers held (not all of them used)
  m_elp = nullptr;
}

Model::~Model()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_random) { delete m_random; m_random = nullptr; }
  if (m_mixModel) { delete m_mixModel; m_mixModel = nullptr; }

#ifndef NDEBUG  // No error checking done and no exceptions thrown in debug mode
  try {
#endif // NDEBUG
    if (m_elp) m_memory->freeEntry(m_elp);
#ifndef NDEBUG
  } catch (...) { cout << "WARNING: Exception in Model::~Model" << endl; }
#endif // NDEBUG
}

void
Model::echo()
//******************************************************************************
//  Echo informaion on model
//! \author  J. Bakosi
//******************************************************************************
{
  // Echo information on selected models
  cout << "Model: " << m_name << "\n";
  cout << " * Number of particles/element: " << m_npel << "\n";
  cout << " * Number of elements: " << m_nel << "\n";
  cout << " * Number of particles: " << m_nel*m_npel << "\n";

  // Echo information on mix model
  cout << " * Mix: "; if (m_mixModel) m_mixModel->echo(); else cout << "none";
  cout << endl;
}

void
Model::init()
//******************************************************************************
//  Initialize model
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_mixModel) m_mixModel->init();
}

void
Model::allocNpel()
//******************************************************************************
//  Allocate array for storing the element IDs of particles
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_elp == nullptr, ModelException,FATAL,ALREADY_ALLOCATED);

  if (m_nel > 1)  // only if inhomogeneous
    m_elp = m_memory->newEntry(m_nel*m_npel, INT, SCALAR, "elp");
}
