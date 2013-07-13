//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.C
  \author    J. Bakosi
  \date      Fri 12 Jul 2013 08:40:14 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Discrete geometry definition
  \details   Discrete geometry definition
*/
//******************************************************************************

#include <DiscreteGeometry.h>
#include <Control.h>

using namespace Quinoa;

DiscreteGeometry::DiscreteGeometry(Memory* const memory,
                                   Paradigm* const paradigm,
                                   Control* const control,
                                   Timer* const timer) :
  Geometry(memory, paradigm, control, timer)
//******************************************************************************
//  Constructor
//! \param[in] memory    Memory oject pointer
//! \param[in] paradigm  Parallel programming paradigm object pointer
//! \param[in] control   Control object
//! \param[in] timer     Timer object
//! \author J. Bakosi
//******************************************************************************
{
  //control->get<control::GEONAME>()
  
}
