//******************************************************************************
/*!
  \file      src/Geometry/AnalyticGeometry.C
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 12:45:16 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Analytic geometry definition
  \details   Analytic geometry definition
*/
//******************************************************************************

#include <AnalyticGeometry.h>
#include <QuinoaControl.h>

using namespace quinoa;

AnalyticGeometry::AnalyticGeometry(const Base& base) : Geometry(base)
//******************************************************************************
//  Constructor
//! \param[in] base  Essentials
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
