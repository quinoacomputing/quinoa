//******************************************************************************
/*!
  \file      src/Model/MixModel/GeneralizedDirichlet/GeneralizedDirichlet.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 01:32:36 PM MST
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

GeneralizedDirichlet::GeneralizedDirichlet(Model* model, const int& nscalar) :
  MixModel(model, "Generalized Dirichlet", nscalar)
//******************************************************************************
//  Constructor
//! \param[in]  model    Model object pointer
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
  cout << "Generalized Dirichlet model:" << endl;
}

void
GeneralizedDirichlet::init()
//******************************************************************************
//  Initialize the generalized Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Initialize generalized Dirichlet model" << endl;
}

void
GeneralizedDirichlet::setIC()
//******************************************************************************
//  Set initial conditions for the generalized Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
}
