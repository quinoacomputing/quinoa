//******************************************************************************
/*!
  \file      src/Model/MixModel/GeneralizedDirichlet/GeneralizedDirichlet.C
  \author    J. Bakosi
  \date      Thu Nov 15 15:53:34 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************

#include <iostream>

#include <GeneralizedDirichlet.h>
#include <MixModel.h>

using namespace std;
using namespace Quinoa;

GeneralizedDirichlet::GeneralizedDirichlet(const int& nscalar) :
  MixModel(nscalar, "Generalized Dirichlet")
//******************************************************************************
//  Constructor
//! \param[in]  nscalar  Number of mixing scalars
//! \author  J. Bakosi
//******************************************************************************
{
}

void
GeneralizedDirichlet::echo()
//******************************************************************************
//  Echo information on the generalized Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  cout << " * Number of mixing scalars: " << m_nscalar << endl;
}

void
GeneralizedDirichlet::init()
//******************************************************************************
//  Initialize the generalized Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  // Set initial conditions
}

void
GeneralizedDirichlet::setIC()
//******************************************************************************
//  Set initial conditions for the generalized Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
}
