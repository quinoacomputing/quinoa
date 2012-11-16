//******************************************************************************
/*!
  \file      src/Model/MixModel/GeneralizedDirichlet/GeneralizedDirichlet.C
  \author    J. Bakosi
  \date      Fri Nov 16 07:33:33 2012
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
