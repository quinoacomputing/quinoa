//******************************************************************************
/*!
  \file      src/Model/Model.C
  \author    J. Bakosi
  \date      Sat 17 Nov 2012 08:09:36 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************

#include <iostream>

#include <Memory.h>
#include <Model.h>

using namespace std;
using namespace Quinoa;

Model::Model(Memory* memory,
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

Model::~Model()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
}
