//******************************************************************************
/*!
  \file      src/Model/Mix/GeneralizedDirichlet/GeneralizedDirichlet.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 01:25:58 PM MST
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

GeneralizedDirichlet::GeneralizedDirichlet(Memory* const memory,
                                           Paradigm* const paradigm,
                                           Control* const control) :
  Mix(memory, paradigm, control, "Generalized Dirichlet")
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
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
