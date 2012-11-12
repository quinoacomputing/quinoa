//******************************************************************************
/*!
  \file      src/Model/MixModel/Dirichlet/Dirichlet.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 01:39:30 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************

#include <iostream>

#include <Dirichlet.h>
#include <MixModel.h>
#include <Model.h>

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
  cout << "Dirichlet model:\n";
  //cout << " * Number of particles/element: " << m_model->npel() << "\n";
  cout << endl;
}

void
Dirichlet::init()
//******************************************************************************
//  Initialize Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Initialize Dirichlet model" << endl;
}

void
Dirichlet::setIC()
//******************************************************************************
//  Set initial conditions for the Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
}
