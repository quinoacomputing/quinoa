//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Sat 19 Jan 2013 11:26:03 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <iostream>

#include <Memory.h>
#include <Physics.h>

using namespace std;
using namespace Quinoa;

Physics::Physics(Memory* memory,
                 Paradigm* paradigm,
                 const string& name,
                 const real time,
                 const int echo,
                 const int nstep) :
  m_memory(memory),
  m_paradigm(paradigm),
  m_name(name),
  m_time(time),
  m_echo(echo),
  m_nstep(nstep)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  name     Name of model
//! \param[in]  time     Maximum time to simulate
//! \param[in]  echo     One-line info in every few time step
//! \param[in]  nstep    Maximum number of time steps to take
//! \author  J. Bakosi
//******************************************************************************
{
}

Physics::~Physics()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
}
