//******************************************************************************
/*!
  \file      src/Model/MixModel/Dirichlet/Dirichlet.C
  \author    J. Bakosi
  \date      Thu Nov 15 15:52:57 2012
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
Dirichlet::init()
//******************************************************************************
//  Initialize Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
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
//   const int npel = m_model->npel();
//   const int nel = m_model->nel();
//   //const VSLStreamStatePtr* stream = m_random->getStr(m_randomStream);
// 
//   int e, p;
// 
//   for (e=0; e<nel; e++ )
//     for (p=0; p<npel; p++ ) {
//       // get a uniformly distributed random numer between [0...1)
//      
//     }
}
