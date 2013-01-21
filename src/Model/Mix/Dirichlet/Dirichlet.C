//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.C
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 11:28:30 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************

#include <iostream>

#include <Dirichlet.h>
#include <Mix.h>

using namespace std;
using namespace Quinoa;

Dirichlet::Dirichlet(const int& nscalar) : Mix(nscalar, "Dirichlet")
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
