//******************************************************************************
/*!
  \file      src/Geometry/Geometry.C
  \author    J. Bakosi
  \date      Mon Oct  7 10:00:43 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry base
  \details   Geometry base
*/
//******************************************************************************

#include <Geometry.h>

using namespace quinoa;

Geometry::Geometry(const Base& base) : m_base(base)
//******************************************************************************
//  Constructor
//! \param[in]  base     Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  //! Echo information on geometry to be created
  echo();
}

void
Geometry::echo()
//******************************************************************************
//  Echo information on geometry
//! \author J. Bakosi
//******************************************************************************
{
  m_base.print.section<ctr::Geometry, ctr::selected, ctr::geometry>();

  m_base.print.subsection("I/O filenames");
  m_base.print.item("Input", m_base.control.get<ctr::io,ctr::input>());
  m_base.print.item("Output", m_base.control.get<ctr::io,ctr::output>());
  m_base.print.endsubsection();
}
