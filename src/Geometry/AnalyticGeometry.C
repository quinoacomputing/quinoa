//******************************************************************************
/*!
  \file      src/Geometry/AnalyticGeometry.C
  \author    J. Bakosi
  \date      Mon Sep  9 08:18:45 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Analytic geometry definition
  \details   Analytic geometry definition
*/
//******************************************************************************

#include <AnalyticGeometry.h>
#include <QuinoaControl.h>

using namespace quinoa;

AnalyticGeometry::AnalyticGeometry(Memory* const memory,
                                   Paradigm* const paradigm,
                                   const QuinoaControl& control,
                                   Timer* const timer) :
  Geometry(memory, paradigm, control, timer),
  m_primitive()
//******************************************************************************
//  Constructor
//! \param[in] memory    Memory oject pointer
//! \param[in] paradigm  Parallel programming paradigm object pointer
//! \param[in] control   Control object
//! \param[in] timer     Timer object
//! \author J. Bakosi
//******************************************************************************
{
  //! Echo information on analytic geometry
  echo();
}

void
AnalyticGeometry::echo()
//******************************************************************************
//  Echo information on analytic geometry
//! \author J. Bakosi
//******************************************************************************
{
  //! Echo information on geometry in general
  Geometry::echo();

  //! Echo information on analytic geometry
}

void
AnalyticGeometry::init()
//******************************************************************************
//  Initialize analytic geometry
//! \author J. Bakosi
//******************************************************************************
{
}

void
AnalyticGeometry::fill()
//******************************************************************************
//  Fill analytic geometry
//! \author J. Bakosi
//******************************************************************************
{
}
