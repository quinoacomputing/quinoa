//******************************************************************************
/*!
  \file      src/Model/MixModel/Dirichlet/Dirichlet.C
  \author    J. Bakosi
  \date      Tue 13 Nov 2012 10:19:11 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************

#include <iostream>

#include <Dirichlet.h>
#include <MixModel.h>
#include <Model.h>
#include <MKLRandom.h>

using namespace std;
using namespace Quinoa;

Dirichlet::Dirichlet(Model* model, const int& nscalar) :
  MixModel(model, "Dirichlet", nscalar)
//******************************************************************************
//  Constructor
//! \param[in]  model    Model object pointer
//! \param[in]  nscalar  Number of mixing scalars
//! \author  J. Bakosi
//******************************************************************************
{
}

void
Dirichlet::echo()
//******************************************************************************
//  Echo information on Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Dirichlet\n";
}

void
Dirichlet::init()
//******************************************************************************
//  Initialize Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  // Allocate data for the Dirichlet model
  m_model->allocNpel();

  // Initialize random number stream
  MKLRandom m_rnd(m_model->memory(), m_model->paradigm());

  // Set initial conditions
  setIC();
}

void
Dirichlet::setIC()
//******************************************************************************
//  Set initial conditions for the Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  const int npel = m_model->npel();
  const int nel = m_model->nel();

  int e, p;

  for (e=0; e<nel; e++ )
    for (p=0; p<npel; p++ ) {
      // get a uniformly distributed random numer between [0...1)
     
    }
}
