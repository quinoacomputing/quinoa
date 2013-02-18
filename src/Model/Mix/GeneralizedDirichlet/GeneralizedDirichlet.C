//******************************************************************************
/*!
  \file      src/Model/Mix/GeneralizedDirichlet/GeneralizedDirichlet.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 10:05:31 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************

#include <iostream>

#include <GeneralizedDirichlet.h>
#include <Mix.h>

using namespace std;
using namespace Quinoa;

GeneralizedDirichlet::GeneralizedDirichlet(Memory* memory,
                                           Paradigm* paradigm,
                                           const int& nscalar,
                                           const int& npar) :
  Mix(memory, paradigm, nscalar, npar, "Generalized Dirichlet")
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  nscalar  Number of mixing scalars
//! \param[in]  npar     Number of particles
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
}
