//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Thu 21 Feb 2013 09:28:58 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <iostream>

//#include <sys/time.h>

#include <Memory.h>
#include <Physics.h>
#include <Control.h>

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
  m_timer(timer)
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
