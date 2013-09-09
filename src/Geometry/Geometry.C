//******************************************************************************
/*!
  \file      src/Geometry/Geometry.C
  \author    J. Bakosi
  \date      Mon Sep  9 08:29:16 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry base
  \details   Geometry base
*/
//******************************************************************************

#include <Geometry.h>
#include <QuinoaControl.h>

using namespace quinoa;

void
Geometry::echo()
//******************************************************************************
//  Echo information on geometry
//! \author J. Bakosi
//******************************************************************************
{
  using namespace control;

  Option<select::Geometry> geo;
  std::cout << " * Geometry: "
            << geo.name(m_control.get<selected,geometry>())
            << std::endl;
}
