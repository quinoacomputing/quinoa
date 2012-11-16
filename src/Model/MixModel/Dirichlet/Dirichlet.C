//******************************************************************************
/*!
  \file      src/Model/MixModel/Dirichlet/Dirichlet.C
  \author    J. Bakosi
  \date      Thu Nov 15 16:29:51 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************

#include <iostream>

#include <Dirichlet.h>
#include <MixModel.h>

using namespace std;
using namespace Quinoa;

Dirichlet::Dirichlet(const int& nscalar) : MixModel(nscalar, "Dirichlet")
//******************************************************************************
//  Constructor
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
  cout << " * Number of mixing scalars: " << m_nscalar << endl;
}

void
Dirichlet::init(const int& npar, real* scalar)
//******************************************************************************
//  Initialize Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  // Set initial conditions
  setIC(npar, scalar);
}

void
Dirichlet::setIC(const int& npar, real* scalar)
//******************************************************************************
//  Set initial conditions for an ensemble of particles
//! \author  J. Bakosi
//******************************************************************************
{
  //const vslstreamstateptr* stream = m_random->getstr(m_randomstream);

  for (int p=0; p<npar; p++ ) {
    // get a uniformly distributed random numer between [0...1)
  }
}
