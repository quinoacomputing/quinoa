//******************************************************************************
/*!
  \file      src/Model/Model.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 01:41:44 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************

#include <iostream>

#include <Model.h>
#include <ModelException.h>
#include <MemoryException.h>
#include <Dirichlet.h>
#include <GeneralizedDirichlet.h>
#include <Setup.h>

using namespace Quinoa;

Model::Model(const ModelType model, const int npel) :
  m_model(model), m_npel(npel)
//******************************************************************************
//  Constructor
//! \param[in]  model  Model type (see Control/Control.h)
//! \param[in]  npel   Number of particles/element
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate model
  // ICC: this could be a switch
  if (model == ModelType::HOMOGENEOUS_DIRICHLET) {
    m_mixModel = new (nothrow) Dirichlet(this,NSCALAR);
    Assert(m_mixModel != nullptr, MemoryException,FATAL,BAD_ALLOC);
  }
  else if (model == ModelType::HOMOGENEOUS_GENDIRICHLET) {
    m_mixModel = new (nothrow) GeneralizedDirichlet(this,NSCALAR);
    Assert(m_mixModel != nullptr, MemoryException,FATAL,BAD_ALLOC);
  } else
    Throw(ModelException,FATAL,NO_SUCH_MODEL);
}

Model::~Model()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_mixModel) { delete m_mixModel; m_mixModel = nullptr; }
}

void
Model::echo()
//******************************************************************************
//  Echo informaion on model
//! \author  J. Bakosi
//******************************************************************************
{
  // Echo information on selected models
  cout << "Model(s):\n";
  cout << " * Mix: ";
  if (m_mixModel) cout << m_mixModel->name(); else cout << "none";
  cout << "\n * Number of particles/element: " << m_npel << "\n";
  cout << endl;

  // Echo information on mix model
  if (m_mixModel) m_mixModel->echo();
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
