//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Wed 06 Mar 2013 06:42:19 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <iostream>

#include <Memory.h>
#include <Paradigm.h>
#include <Control.h>
#include <Physics.h>

using namespace std;
using namespace Quinoa;
using namespace control;

Physics::Physics(Memory* const memory,
                 Paradigm* const paradigm,
                 Control* const control,
                 Timer* const timer) :
  m_memory(memory),
  m_paradigm(paradigm),
  m_control(control),
  m_timer(timer),
  m_nthread(paradigm->nthread())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object
//! \param[in]  paradigm Parallel programming object
//! \param[in]  control  Control object
//! \param[in]  timer    Timer object
//! \author  J. Bakosi
//******************************************************************************
{
}
